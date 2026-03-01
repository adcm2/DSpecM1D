#%%
import numpy as np
import matplotlib.pyplot as plt
import math
import re


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
fig, axes = plt.subplots(1, 1, figsize=(14, 12), sharex=True)
plt.style.use('seaborn-v0_8-whitegrid')

# --- 3. Define Plotting Parameters ---
lwidth = 2.0
BIGGER_SIZE = 20
# time_vector = dmf[:, 0]

# find lower and upper indices for frequency range of interest (0.1 to 6 mHz)
fmin, fmax = 0.0, 5.0
idx_min = np.searchsorted(dmf[:, 0], fmin)
idx_max = np.searchsorted(dmf[:, 0], fmax, side='right')
time_vector = dmf[idx_min:idx_max, 0]
print(f"Plotting frequency range: {dmf[idx_min, 0]:.2f} mHz to {dmf[idx_max,0]:.2f} mHz (indices {idx_min} to {idx_max-1})")
# --- Data Column Indices ---
# Give meaningful names to the columns being plotted.
# This makes the code MUCH easier to read and debug.
trans_z, trans_n, trans_e = dmf[idx_min:idx_max, 3], dmf[idx_min:idx_max, 6], dmf[idx_min:idx_max, 9]
yspec_z, yspec_n, yspec_e = dmf[idx_min:idx_max, 12], dmf[idx_min:idx_max, 15], dmf[idx_min:idx_max, 18]
mineos_z, mineos_n, mineos_e = dmf[idx_min:idx_max, 21], dmf[idx_min:idx_max, 24], dmf[idx_min:idx_max, 27]


# --- 4. Calculate Normalization Factors ---
# Find the peak absolute amplitude for each YSpec component to use for normalization.
# Add a small epsilon to prevent division by zero if a trace is flat.
norm_z = np.max(np.abs(yspec_z)) + 1e-12
norm_n = np.max(np.abs(yspec_n)) + 1e-12
norm_e = np.max(np.abs(yspec_e)) + 1e-12

# --- 5. Plot Data and Relative Differences for Each Component ---

# == Z Component (Row 0) ==
ax_data = axes
ax_data.plot(time_vector, yspec_z, "b", linewidth=lwidth, label='YSpec')
ax_data.plot(time_vector, trans_z, "r--", linewidth=lwidth, label='DSpecM1D')
# ax_data.plot(time_vector, mineos_z, "g-.", linewidth=lwidth, label='MINEOS')

# set labels and legend
ax_data.set_xlabel('Frequency (mHz)', fontsize=BIGGER_SIZE)
ax_data.set_ylabel('Z displacement (m)', fontsize=BIGGER_SIZE)
ax_data.legend(loc='upper left', fontsize = BIGGER_SIZE)
ax_data.tick_params(axis='both', which='major', labelsize=BIGGER_SIZE)
ax_data.set_xlim([fmin, fmax])

# calculate relative difference
yspec_diff = np.abs(yspec_z - trans_z)/norm_z * 100
yspec_av_diff = np.mean(yspec_diff)
# mineos_diff = np.abs(mineos_z - trans_z)/norm_z * 100
# mineos_av_diff = np.mean(mineos_diff)

print(f"Average relative difference between YSpec and Transverse (Z): {yspec_av_diff:.4f} %")
# print(f"Average relative difference between MINEOS and Transverse (Z): {mineos_av_diff:.4f} %")

# print maximum relative differences for more insight
print(f"Maximum relative difference between YSpec and Transverse (Z): {np.max(yspec_diff):.4f} %")
# print(f"Maximum relative difference between MINEOS and Transverse (Z): {np.max(mineos_diff):.4f} %")


# Improve layout
plt.tight_layout(rect=[0, 0, 1, 0.94]) # Adjust rect to make space for suptitle
plt.show()