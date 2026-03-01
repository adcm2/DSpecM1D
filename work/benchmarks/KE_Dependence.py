#%%
import numpy as np
import matplotlib.pyplot as plt
import math
import re

# This helper function is fine, no changes needed here.
_float_re = re.compile(r'[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][+-]?\d+)?')
def read_spectra_to_array(path, pad=np.nan):
    # ... (existing function is ok)
    rows = []
    with open(path, 'r') as fh:
        for line in fh:
            nums = _float_re.findall(line)
            if not nums:
                continue
            rows.append([float(x) for x in nums])
    if not rows:
        return np.empty((0, 0), dtype=float)
    maxcols = max(len(r) for r in rows)
    arr = np.full((len(rows), maxcols), pad, dtype=float)
    for i, r in enumerate(rows):
        arr[i, :len(r)] = r
    return arr

# --- 1. Load Data ---
# Use a try-except block for better error handling.
try:
    path_freqs = "tidal_frequency.out"
    dfreqs = np.loadtxt(path_freqs, delimiter=";")
except FileNotFoundError:
    print(f"Error: Data file not found at '{path_mf}'")
    exit()

# --- 2. Setup a Cleaner Figure ---
# Create a 3x2 grid. Left column for data, right column for differences.
# `sharex=True` links the x-axes for easier comparison.
fig, axes = plt.subplots(1, 1, figsize=(14, 12))
plt.style.use('seaborn-v0_8-whitegrid')

# --- 3. Define Plotting Parameters ---
lw = 1.5
BIGGER_SIZE = 14

# --- Data Column Indices ---
# Give meaningful names to the columns being plotted.
# This makes the code MUCH easier to read and debug.
per_val = np.log(dfreqs[:, 0])
ke_ratio = np.log(dfreqs[:, 1])
plt.plot(per_val, ke_ratio, "b", linewidth=lw)
plt.title('Kinetic Energy Ratio vs Period', fontsize=BIGGER_SIZE)
plt.xlabel('Period (BV units)', fontsize=BIGGER_SIZE)


# Improve layout
# plt.tight_layout(rect=[0, 0, 1, 0.94]) # Adjust rect to make space for suptitle
plt.show()