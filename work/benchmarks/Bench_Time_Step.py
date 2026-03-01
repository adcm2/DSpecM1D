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
    path_mf = "bench_t_step.out"
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
nz = int((dmf.shape[1] - 1)/3)
sd = 3 * (nz - 1)
z_1, n_1, e_1 = dmf[:, 1], dmf[:, 2], dmf[:, 3]
z_2, n_2, e_2 = dmf[:, 4], dmf[:, 5] , dmf[:, 6]   # Adjust East component
z_3, n_3, e_3 = dmf[:, 7], dmf[:, 8], dmf[:, 9]
z_4, n_4, e_4 = dmf[:, 1 + sd], dmf[:, 2 + sd], dmf[:, 3 + sd]


# --- 4. Calculate Normalization Factors ---
# Find the peak absolute amplitude for each YSpec component to use for normalization.
# Add a small epsilon to prevent division by zero if a trace is flat.
norm_z = np.max(np.abs(z_2)) + 1e-12
norm_n = np.max(np.abs(n_2)) + 1e-12
norm_e = np.max(np.abs(e_2)) + 1e-12

# --- 5. Plot Data and Relative Differences for Each Component ---

xlow = 400
xtop = 1300

# == Z Component (Row 0) ==
ax_data = axes[0, 0]
ax_diff = axes[0, 1]

for i in range(0,nz,4):
    ax_data.plot(time_vector, dmf[:, 1 + 3 * i],  linewidth=lwidth, label=f'YSpec {i+1}')
ax_data.plot(time_vector, z_4, "k", linewidth=lwidth, label='0.001')
# ax_data.plot(time_vector, z_3, "b", linewidth=lwidth, label='0.01')
# ax_data.plot(time_vector, z_2, "r--", linewidth=lwidth, label='0.05')
# ax_data.plot(time_vector, z_1, "g-.", linewidth=lwidth, label='0.1')
ax_data.set_ylabel('Z-Component', fontsize=BIGGER_SIZE)
ax_data.legend(loc='upper right')
# ax_data.set_xlim([xlow,xtop])

# Calculate and plot relative difference
ax_diff.plot(time_vector, (z_1 - z_4) / norm_z * 100, "k", linewidth=lwidth, label='Error 0.1')
ax_diff.plot(time_vector, (z_2 - z_4) / norm_z * 100, "m", linewidth=lwidth, label='Error 0.05')
ax_diff.plot(time_vector, (z_3 - z_4) / norm_z * 100, "b", linewidth=lwidth, label='Error 0.01')
ax_diff.set_ylabel('Relative Error (%)', fontsize=BIGGER_SIZE)
ax_diff.legend(loc='upper right')
# ax_diff.set_xlim([xlow,xtop])

# == North Component (Row 1) ==
ax_data = axes[1, 0]
ax_diff = axes[1, 1]

ax_data.plot(time_vector, n_4, "k", linewidth=lwidth)
ax_data.plot(time_vector, n_3 , "b", linewidth=lwidth)
ax_data.plot(time_vector, n_2, "r--", linewidth=lwidth)
ax_data.plot(time_vector, n_1, "g-.", linewidth=lwidth)
ax_data.set_ylabel('North-Component', fontsize=BIGGER_SIZE)
# ax_data.set_xlim([xlow,xtop])

# Calculate and plot relative difference
ax_diff.plot(time_vector, (n_1 - n_4) / norm_n * 100, "k", linewidth=lwidth, label='Error 0.1')
ax_diff.plot(time_vector, (n_2 - n_4) / norm_n * 100, "m", linewidth=lwidth, label='Error 0.05')
ax_diff.plot(time_vector, (n_3 - n_4) / norm_n * 100, "b", linewidth=lwidth, label='Error 0.01')
ax_diff.set_ylabel('Relative Error (%)', fontsize=BIGGER_SIZE)
# ax_diff.set_xlim([xlow,xtop])

# == East Component (Row 2) ==
ax_data = axes[2, 0]
ax_diff = axes[2, 1]

ax_data.plot(time_vector, e_4, "k", linewidth=lwidth)
ax_data.plot(time_vector, e_3, "b", linewidth=lwidth)
ax_data.plot(time_vector, e_2, "r--", linewidth=lwidth);
ax_data.plot(time_vector, e_1, "g-.", linewidth=lwidth);

ax_data.set_ylabel('East-Component', fontsize=BIGGER_SIZE)
# ax_data.set_xlim([xlow,xtop])

# Calculate and plot relative difference
ax_diff.plot(time_vector, (e_1 - e_4) / norm_e * 100, "k", linewidth=lwidth)
ax_diff.plot(time_vector, (e_2 - e_4) / norm_e * 100, "m", linewidth=lwidth)
ax_diff.plot(time_vector, (e_3 - e_4) / norm_e * 100, "b", linewidth=lwidth)
# ax_diff.set_xlim([xlow,xtop])
ax_diff.set_ylabel('Relative Error (%)', fontsize=BIGGER_SIZE)

# --- 6. Final Figure Adjustments ---
# Set titles for the columns
axes[0, 0].set_title('Seismic Traces', fontsize=BIGGER_SIZE + 2)
axes[0, 1].set_title('Relative Difference (%)', fontsize=BIGGER_SIZE + 2)

# Set shared x-axis properties for the bottom plots
axes[2, 0].set_xlabel('Time (s)', fontsize=BIGGER_SIZE)
axes[2, 1].set_xlabel('Time (s)', fontsize=BIGGER_SIZE)
# axes[2, 0].set_xlim(400, 1250)

# Add a main title to the entire figure
fig.suptitle('Comparison of YSpec and Transverse Models (High Frequency)', fontsize=BIGGER_SIZE + 4, y=0.96)

# Improve layout
plt.tight_layout(rect=[0, 0, 1, 0.94]) # Adjust rect to make space for suptitle
plt.show()