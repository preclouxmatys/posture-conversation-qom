"""
01_compute_qom_vicon.py

Computes the Quantity of Movement (QoM) of the head and wrists for each
participant x posture recording, from the raw Vicon CSV exports.

Input : data/raw/vicon_csv/<CONDITION><DYAD>.csv   (e.g. SEATEDD01.csv)
        CONDITION in {SEATED, SEMI, STANDING}; one file per dyad and posture,
        containing the markers of both participants (P2 labels are prefixed "2").
        These files are not shared (see README, "Data availability").
Output: data/final/posture_qom_long.xlsx
        long format, 120 rows (20 dyads x 3 postures x 2 participants)

Processing, per participant and recording (see Section 2.6 of the article):
  1. Brief marker-label swaps are detected from aberrant segment lengths
     (shoulder-elbow or elbow-wrist length deviating by more than 25% from its
     median); the distal marker is set to missing over the affected samples
     (+/- 5 samples).
  2. Gaps of up to 2 s are filled by linear interpolation.
  3. Trajectories are low-pass filtered (4th-order Butterworth, 6 Hz, applied
     forward and backward: zero phase lag).
  4. Head position = midpoint of the two temple markers, computed at each sample
     (filtered markers), so that head QoM is the path length of a single point.
  5. QoM = cumulative 3D Euclidean displacement between successive samples,
     summed over both wrists for the upper limbs.
  6. Normalization by the participant's median shoulder width in the recording
     (distance between the two shoulder markers): QoM is expressed in shoulder widths.
Positions are expressed in the laboratory reference frame.
"""
import csv
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt

import config

CONDITIONS = ["SEATED", "SEMI", "STANDING"]
MARKERS = ["Tempe_D", "Tempe_G", "epaule_D", "epaule_G",
           "coude_D", "coude_G", "poignet_D", "poignet_G"]
P2_ALIAS = {"Tempe_G": "Temps_G"}  # label typo for participant 2 in some raw files
OUTPUT_PATH = Path(config.FINAL_DIR) / "posture_qom_long.xlsx"


# ------------------------------------------------------------------ reading
def read_vicon(path):
    """Reads the two Vicon CSV export formats present in the data
    (header with 'Frame, Sub Frame' and a units line, or 'Field #' without
    Sub Frame or units)."""
    with open(path, encoding="utf-8", errors="replace", newline="") as f:
        rd = csv.reader(f)
        rows = [next(rd) for _ in range(6)]
    ax_i = next(i for i, r in enumerate(rows) if r and r[0].strip() in ("Frame", "Field #"))
    mk, ax = rows[ax_i - 1], rows[ax_i]
    skip = ax_i + 1 + (1 if any(c.strip() == "mm" for c in rows[ax_i + 1]) else 0)
    names, last = [], ""
    mk = mk + [""] * (len(ax) - len(mk))
    for m, a in zip(mk, ax):
        last = m.strip() or last
        a = a.strip()
        if a in ("X", "Y", "Z"):
            names.append(f"{last.split(':')[-1]}_{a}")
        else:
            names.append("Frame" if a in ("Frame", "Field #") else (a or f"_c{len(names)}"))
    df = pd.read_csv(path, skiprows=skip, header=None, engine="c", low_memory=False)
    df = df.iloc[:, :len(names)]
    df.columns = names[:df.shape[1]]
    df = df.apply(pd.to_numeric, errors="coerce")
    return df[df["Frame"].notna()].reset_index(drop=True)


def get_marker(df, name):
    cols = [f"{name}_{a}" for a in "XYZ"]
    if not all(c in df.columns for c in cols):
        return None
    return df[cols].to_numpy(float)


def participant_markers(df, pid):
    prefix = "" if pid == "P1" else "2"
    out = {}
    for m in MARKERS:
        arr = get_marker(df, prefix + (P2_ALIAS.get(m, m) if pid == "P2" else m))
        if arr is None and pid == "P2":
            arr = get_marker(df, prefix + m)
        out[m] = arr if arr is not None else np.full((len(df), 3), np.nan)
    return out


# ------------------------------------------------------------------ processing
def interp_nan(x, maxgap=None):
    """Linear interpolation of missing samples; gaps longer than maxgap stay missing."""
    df = pd.DataFrame(x)
    if maxgap is None:
        return df.interpolate(limit_direction="both").to_numpy(copy=True)
    out = df.interpolate(limit_direction="both", limit=maxgap, limit_area=None).to_numpy(copy=True)
    isn = np.isnan(np.asarray(x, float)).any(axis=1) if np.ndim(x) > 1 else np.isnan(x)
    run = np.zeros(len(isn), int)
    c = 0
    for i, v in enumerate(isn):
        c = c + 1 if v else 0
        run[i] = c
    longgap = np.zeros(len(isn), bool)
    i = len(isn) - 1
    while i >= 0:
        if run[i] > maxgap:
            L = run[i]
            longgap[i - L + 1:i + 1] = True
            i -= L
        else:
            i -= 1
    out[longgap] = np.nan
    return out


def lowpass(x, fs):
    """Zero-phase Butterworth low-pass filter; samples missing before filtering stay missing."""
    x = np.asarray(x, float)
    mask = np.isnan(x).any(axis=1)
    b, a = butter(config.BUTTER_ORDER, config.CUTOFF_HZ / (fs / 2))
    y = filtfilt(b, a, interp_nan(x), axis=0)
    y[mask] = np.nan
    return y


def remove_marker_swaps(V):
    """Brief label swaps: shoulder-elbow or elbow-wrist length deviating by more
    than SWAP_TOLERANCE from its median -> distal marker set to missing."""
    n_bad = 0
    width = 2 * config.SWAP_DILATION_SAMPLES + 1
    for s in "DG":
        for prox, dist in [(f"epaule_{s}", f"coude_{s}"), (f"coude_{s}", f"poignet_{s}")]:
            L = np.linalg.norm(V[dist] - V[prox], axis=1)
            med = np.nanmedian(L)
            bad = np.abs(L - med) > config.SWAP_TOLERANCE * med
            bad = np.convolve(bad.astype(int), np.ones(width), "same") > 0
            V[dist][bad] = np.nan
            n_bad += int(bad.sum())
    return n_bad


def path_length(p):
    return float(np.nansum(np.linalg.norm(np.diff(p, axis=0), axis=1)))


def process(V):
    n_swap = remove_marker_swaps(V)
    Vf = {}
    for k, a in V.items():
        if np.isfinite(a).mean() < 0.5:  # marker missing for most of the recording
            Vf[k] = np.full_like(a, np.nan)
            continue
        Vf[k] = lowpass(interp_nan(a, maxgap=config.MAX_GAP_SAMPLES), config.FS_VICON)
    head = (Vf["Tempe_D"] + Vf["Tempe_G"]) / 2
    shoulder_width = float(np.nanmedian(np.linalg.norm(Vf["epaule_D"] - Vf["epaule_G"], axis=1)))
    qom_wrist_mm = path_length(Vf["poignet_D"]) + path_length(Vf["poignet_G"])
    qom_head_mm = path_length(head)
    return dict(shoulder_width_mm=shoulder_width,
                QoM_wrist_mm=qom_wrist_mm, QoM_head_mm=qom_head_mm,
                QoM_wrist_norm=qom_wrist_mm / shoulder_width,
                QoM_head_norm=qom_head_mm / shoulder_width,
                n_samples=len(head), n_swap_samples=n_swap)


def main():
    files = sorted(Path(config.VICON_CSV_DIR).glob("*.csv"))
    if not files:
        raise FileNotFoundError(f"No Vicon CSV file in {config.VICON_CSV_DIR} (see README).")
    rows = []
    for path in files:
        cond = next(c for c in CONDITIONS if path.stem.upper().startswith(c))
        dyad = path.stem[len(cond):]
        df = read_vicon(path)
        for pid in ["P1", "P2"]:
            res = process(participant_markers(df, pid))
            rows.append(dict(dyad=dyad, condition=cond, pid=pid, **res))
            print(f"{dyad} {cond} {pid}: {res['n_samples']} samples, "
                  f"wrist QoM = {res['QoM_wrist_norm']:.2f}, head QoM = {res['QoM_head_norm']:.2f}")

    out = pd.DataFrame(rows).sort_values(["dyad", "condition", "pid"]).reset_index(drop=True)
    out["posture"] = pd.Categorical(out["condition"], categories=CONDITIONS, ordered=True)
    out["subject"] = out["dyad"] + "_" + out["pid"]
    assert out["dyad"].nunique() == 20 and out["subject"].nunique() == 40 and len(out) == 120, \
        "Expected 20 dyads, 40 participants and 120 recordings"
    cols = ["dyad", "condition", "pid", "shoulder_width_mm", "QoM_wrist_mm", "QoM_head_mm",
            "QoM_wrist_norm", "QoM_head_norm", "posture", "subject"]
    out[cols].to_excel(OUTPUT_PATH, index=False)
    print(f"Saved: {OUTPUT_PATH} ({len(out)} rows)")


if __name__ == "__main__":
    main()
