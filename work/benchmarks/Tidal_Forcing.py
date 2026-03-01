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
    path_mf = "tidal_radial_response.out"
    dmf = np.loadtxt(path_mf, delimiter=";")
    path_freqs = "tidal_frequency.out"
    dfreqs = np.loadtxt(path_freqs, delimiter=";")
    path_n2 = "tidal_N2.out"
    dn2 = np.loadtxt(path_n2, delimiter=";")
except FileNotFoundError:
    print(f"Error: Data file not found at '{path_mf}'")
    exit()

# --- 2. Setup a Cleaner Figure ---
# Create a 3x2 grid. Left column for data, right column for differences.
# `sharex=True` links the x-axes for easier comparison.
nfreqs = dfreqs.shape[0]  # Number of frequency columns
print(f"Number of frequencies found: {nfreqs}")
if(nfreqs > 6):
    nfreqs = 6  # Limit to first 6 frequencies for plotting
fig, axes = plt.subplots(1, nfreqs + 1, figsize=(14, 12))
plt.style.use('seaborn-v0_8-whitegrid')

# --- 3. Define Plotting Parameters ---
lw = 1.5
BIGGER_SIZE = 14
time_vector = dmf[:, 0]

# --- Data Column Indices ---
# Give meaningful names to the columns being plotted.
# This makes the code MUCH easier to read and debug.
rad_val = dmf[:, 0]
u_real = dmf[:, 1]
u_imag = dmf[:, 2]
u_abs = dmf[:, 3]
v_real = dmf[:, 4]
v_imag = dmf[:, 5]
v_abs = dmf[:, 6]
# u_abs = dmf[:, 1]
# v_abs = dmf[:, 2]
low_rad = dn2[0,0]
up_rad = dn2[-1,0]
for i in range(nfreqs):
    freq_hours = dfreqs[i]
    ax_data = axes[i]
    ax_data.plot(dmf[:,6 * i + 1], rad_val,"b", linewidth=lw, label='U real')
    ax_data.plot(dmf[:,6 * i + 2], rad_val, "r", linewidth=lw, label='U imag')
    ax_data.plot(dmf[:,6 * i + 4], rad_val,"g", linewidth=lw, label='V real')
    ax_data.plot(dmf[:,6 * i + 5], rad_val, "m", linewidth=lw, label='V imag')
    ax_data.set_title(dfreqs[i,0], fontsize=BIGGER_SIZE)
    ax_data.set_ylabel('Radius (m)', fontsize=BIGGER_SIZE)
    ax_data.set_xlabel('Amplitude', fontsize=BIGGER_SIZE)
    ax_data.tick_params(axis='both', which='major', labelsize=BIGGER_SIZE - 2)
    ax_data.legend(loc='upper right')
    ax_data.set_ylim([low_rad, up_rad])

ax_data = axes[nfreqs]
ax_data.plot(dn2[:,1], dn2[:,0],"g", linewidth=lw)
ax_data.axvline(0, color='k', linestyle='--', linewidth=1)
ax_data.set_title('N² Profile', fontsize=BIGGER_SIZE)
ax_data.set_ylim([low_rad, up_rad])

# Improve layout
plt.tight_layout(rect=[0, 0, 1, 0.94]) # Adjust rect to make space for suptitle
plt.show()