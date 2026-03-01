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
    path_mf = "bench_w_lf.out"
    dmf = np.loadtxt(path_mf, delimiter=";")
except FileNotFoundError:
    print(f"Error: Data file not found at '{path_mf}'")
    exit()

# --- 2. Setup a Cleaner Figure ---
# Create a 3x2 grid. Left column for data, right column for differences.
# `sharex=True` links the x-axes for easier comparison.
fig, axes = plt.subplots(3, 2, figsize=(14, 12), sharex=True)
plt.style.use('seaborn-v0_8-whitegrid')

# --- 3. Define Plotting Parameters ---
lwidth = 1.5
BIGGER_SIZE = 14
time_vector = dmf[:, 0]

# --- Data Column Indices ---
# Give meaningful names to the columns being plotted.
# This makes the code MUCH easier to read and debug.
trans_z, trans_n, trans_e = dmf[:, 3], dmf[:, 6], dmf[:, 9]
yspec_z, yspec_n, yspec_e = dmf[:, 12], dmf[:, 15], dmf[:, 18]
mineos_z, mineos_n, mineos_e = dmf[:, 21], dmf[:, 24], dmf[:, 27]


# --- 4. Calculate Normalization Factors ---
# Find the peak absolute amplitude for each YSpec component to use for normalization.
# Add a small epsilon to prevent division by zero if a trace is flat.
norm_z = np.max(np.abs(yspec_z)) + 1e-12
norm_n = np.max(np.abs(yspec_n)) + 1e-12
norm_e = np.max(np.abs(yspec_e)) + 1e-12

# --- 5. Plot Data and Relative Differences for Each Component ---

# == Z Component (Row 0) ==
ax_data = axes[0, 0]
ax_diff = axes[0, 1]

ax_data.plot(time_vector, yspec_z, "b", linewidth=lwidth, label='YSpec')
ax_data.plot(time_vector, trans_z, "r--", linewidth=lwidth, label='Transverse')
# ax_data.plot(time_vector, mineos_z, "g-.", linewidth=lwidth, label='MINEOS')
ax_data.set_ylabel('Z-Component', fontsize=BIGGER_SIZE)
ax_data.legend(loc='upper right')

# Calculate and plot relative difference
ax_diff.plot(time_vector, (trans_z - yspec_z) / norm_z * 100, "k", linewidth=lwidth, label='(Trans - YSpec) / Peak(YSpec)')
# ax_diff.plot(time_vector, (trans_z - mineos_z) / norm_z * 100, "purple", linestyle='--', linewidth=lwidth, label='(Trans - MINEOS) / Peak(YSpec)')
ax_diff.set_ylabel('Relative Error (%)', fontsize=BIGGER_SIZE)
ax_diff.legend(loc='upper right')

# == North Component (Row 1) ==
ax_data = axes[1, 0]
ax_diff = axes[1, 1]

ax_data.plot(time_vector, dmf[:,6], "b", linewidth=lwidth)
ax_data.plot(time_vector, dmf[:,15], "r--", linewidth=lwidth)
# ax_data.plot(time_vector, mineos_n, "g-.", linewidth=lwidth)
ax_data.set_ylabel('North-Component', fontsize=BIGGER_SIZE)

# Calculate and plot relative difference
ax_diff.plot(time_vector, (trans_n - yspec_n) / norm_n * 100, "k", linewidth=lwidth)
# ax_diff.plot(time_vector, (trans_n - mineos_n) / norm_n * 100, "purple", linestyle='--', linewidth=lwidth)
ax_diff.set_ylabel('Relative Error (%)', fontsize=BIGGER_SIZE)

# == East Component (Row 2) ==
ax_data = axes[2, 0]
ax_diff = axes[2, 1]

ax_data.plot(time_vector, dmf[:,9], "b", linewidth=lwidth)
ax_data.plot(time_vector, dmf[:,18], "r--", linewidth=lwidth)
# ax_data.plot(time_vector, mineos_e, "g-.", linewidth=lwidth)
ax_data.set_ylabel('East-Component', fontsize=BIGGER_SIZE)

# Calculate and plot relative difference
ax_diff.plot(time_vector, (trans_e - yspec_e) / norm_e * 100, "k", linewidth=lwidth)
# ax_diff.plot(time_vector, (trans_e - mineos_e) / norm_e * 100, "purple", linestyle='--', linewidth=lwidth)
ax_diff.set_ylabel('Relative Error (%)', fontsize=BIGGER_SIZE)

# --- 6. Final Figure Adjustments ---
# Set titles for the columns
axes[0, 0].set_title('Seismic Traces', fontsize=BIGGER_SIZE + 2)
axes[0, 1].set_title('Relative Difference (%)', fontsize=BIGGER_SIZE + 2)

# Set shared x-axis properties for the bottom plots
axes[2, 0].set_xlabel('Frequency (mHz)', fontsize=BIGGER_SIZE)
axes[2, 1].set_xlabel('Frequency (mHz)', fontsize=BIGGER_SIZE)
# axes[2, 0].set_xlim(400, 1250)

# Add a main title to the entire figure
fig.suptitle('Comparison of YSpec and Transverse Models (High Frequency)', fontsize=BIGGER_SIZE + 4, y=0.96)

# Improve layout
plt.tight_layout(rect=[0, 0, 1, 0.94]) # Adjust rect to make space for suptitle
plt.show()