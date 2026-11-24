"""Paths and processing parameters (Python side; see config_r.R for the R side)."""
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Raw Vicon exports (identifiable data, never versioned; see README, "Data availability")
VICON_CSV_DIR = os.path.join(BASE_DIR, "data", "raw", "vicon_csv")

# Versioned outputs
FINAL_DIR = os.path.join(BASE_DIR, "data", "final")
FIG_DIR = os.path.join(FINAL_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Acquisition
FS_VICON = 100.0  # Hz

# Gap filling: linear interpolation of gaps up to this many samples (2 s at 100 Hz)
MAX_GAP_SAMPLES = 200

# Low-pass filter: Butterworth, zero-phase (filtfilt)
BUTTER_ORDER = 4
CUTOFF_HZ = 6.0

# Marker-swap detection: segment length deviating from its median by more than this fraction
SWAP_TOLERANCE = 0.25
SWAP_DILATION_SAMPLES = 5  # samples flagged on each side of a detected swap

# Reproducibility (bootstrap)
RANDOM_SEED = 42
N_BOOT = 5000
