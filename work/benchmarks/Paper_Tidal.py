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
    path_mf = "tidal_radial_response0.001000.out"
    dmf = np.loadtxt(path_mf, delimiter=";")
    
    path_freqs = "tidal_frequency_0.001000.out"
    dfreqs = np.loadtxt(path_freqs, delimiter=";")
    
    path_n2 = "tidal_N2_0.001000.out"
    dn2 = np.loadtxt(path_n2, delimiter=";")
    
except FileNotFoundError:
    print(f"Error: Data file not found at '{path_mf}'")
    exit()

# --- 2. Setup a Cleaner Figure ---
# Create a 1-row grid.
nfreqs = dfreqs.shape[0]  # Number of frequency columns
print(f"Number of frequencies found: {nfreqs}")
if(nfreqs > 6):
    nfreqs = 6  # Limit to first 6 frequencies for plotting
    
# Changed to 1 row
fig, axes = plt.subplots(1, nfreqs + 1, figsize=(14, 6), sharey=True)
fig.subplots_adjust(wspace=0) # Remove horizontal space between plots
plt.style.use('seaborn-v0_8-whitegrid')

# --- 3. Define Plotting Parameters ---
lw = 1.5
M_SIZE = 14
BIGGER_SIZE = 18
time_vector = dmf[:, 0]

# --- Data Column Indices ---
# Give meaningful names to the columns being plotted.
# This makes the code MUCH easier to read and debug.
rad_val = dmf[:, 0]
# rad_val2 removed
u_abs = dmf[:, 1]
v_abs = dmf[:, 2]
low_rad = dn2[0,0]
up_rad = dn2[-1,0]
# List of Roman numerals for labeling
roman_numerals = ['i)', 'ii)', 'iii)', 'iv)', 'v)', 'vi)']

# Loop for the single row
for i in range(nfreqs):
    freq_hours = dfreqs[i]
    # Access axis directly since it's now 1D array
    ax_data = axes[i]
    maxval_u = np.max(np.abs(dmf[:, 2 * i + 1]))
    maxval_v = np.max(np.abs(dmf[:, 2 * i + 2]))
    # print(f"Frequency: {freq_hours} hours, Max U: {maxval_u}, Max V: {maxval_v}")#
    maxvalt = 1.1*max(np.max(np.abs(dmf[:, 6 * i + 1])), np.max(np.abs(dmf[:, 6 * i + 2])),np.max(np.abs(dmf[:, 6 * i + 4])), np.max(np.abs(dmf[:, 6 * i + 5])))
    ax_data.plot(dmf[:,6 * i + 1], rad_val,"b", linewidth=lw, label='U')
    # ax_data.plot(dmf[:,6 * i + 2], rad_val, "r", linewidth=lw, label='U imag')
    ax_data.plot(dmf[:,6 * i + 4], rad_val,"r", linewidth=lw, label='V')
    # ax_data.plot(dmf[:,6 * i + 5], rad_val, "m", linewidth=lw, label='V imag')
    
    # Add text label i) freq, ii) freq etc. in top left
    # x=0.05, y=0.95 puts it in the top-left corner, relative to axis size
    label_text = f"{roman_numerals[i]} {dfreqs[i,0]:.1f}"
    ax_data.text(0.05, 0.95, label_text, transform=ax_data.transAxes, 
                 fontsize=M_SIZE, verticalalignment='top', fontweight='bold')

    # Add horizontal lines at limits
    ax_data.axhline(low_rad, color='k', linestyle='--', linewidth=1)
    ax_data.axhline(up_rad, color='k', linestyle='--', linewidth=1)
    
    # Add vertical line at x=0
    ax_data.axvline(0, color='k', linestyle='--', linewidth=1)

    # Thicken x-axis
    # ax_data.spines['bottom'].set_linewidth(2.0)
    ax_data.set_ylim([0, rad_val[-1]]) # Set y-limits to match the data range
    ax_data.set_xlim([-maxvalt, maxvalt]) # Set x-limits to be symmetric around zero based on max value
    # ax_data.set_ylabel('Radius (m)', fontsize=BIGGER_SIZE)
    # ax_data.set_xlabel('Amplitude', fontsize=BIGGER_SIZE)
    # ax_data.tick_params(axis='both', which='major', labelsize=BIGGER_SIZE - 2)
    ax_data.set_xticks([0]) 
    ax_data.set_xticklabels([0]) # Remove x-axis labels for all but the first plot
    ax_data.spines['left'].set_linewidth(1.5)
    ax_data.spines['right'].set_linewidth(1.5)
    ax_data.spines['bottom'].set_linewidth(1.5)
    ax_data.spines['top'].set_linewidth(1.5)
    ax_data.tick_params(axis='x', which='major', length=10, width=1.5)
    ax_data.tick_params(axis='x', which='major', labelcolor='black',labelsize=M_SIZE)
   # Remove x-axis ticks for all but the first plot
    if i > 0: # Remove y-axis labels for all but the first plot
        ax_data.set_yticklabels([])
        ax_data.set_yticks([])
        
        # ax_data.spines['left'].set_visible(False)
        # ax_data.spines['right'].set_visible(False)

ax_data=axes[0]
# Bundle weight and size into the properties dictionary for consistent behavior
legend_properties = {'weight': 'bold', 'size': BIGGER_SIZE}

# ax_data.legend(loc='upper center', bbox_to_anchor=(2.5, 1.1), 
            #    ncol=4, fancybox=True, shadow=True, 
            #    frameon=True, edgecolor='black', prop=legend_properties)

# N^2 Plot (now at the end of the single row)
ax_data = axes[nfreqs]
ax_data.plot(dn2[:,1], dn2[:,0],"g", linewidth=lw)
ax_data.axvline(0, color='k', linestyle='--', linewidth=1)
ax_data.spines['left'].set_linewidth(1.5)
ax_data.spines['right'].set_linewidth(1.5)
ax_data.spines['bottom'].set_linewidth(1.5)
ax_data.spines['top'].set_linewidth(1.5)
ax_data.tick_params(axis='x', which='major', length=10, width=1.5)
ax_data.tick_params(axis='x', which='major', labelcolor='black',labelsize=M_SIZE)
# Add horizontal lines at limits
ax_data.axhline(low_rad, color='k', linestyle='--', linewidth=1)
ax_data.axhline(up_rad, color='k', linestyle='--', linewidth=1)
ax_data.set_xticks([0]) 
ax_data.set_xticklabels([0]) # Remove x-axis labels for all but the first 
ax_data.text(0.05, 0.95, "v) N² Profile", transform=ax_data.transAxes, 
                 fontsize=M_SIZE, verticalalignment='top', fontweight='bold')
ax_data.set_yticklabels([])
ax_data.set_yticks([])

# Improve layout
# Call tight_layout FIRST to fit the outer labels, then zero out the wspace
# plt.tight_layout(rect=[0, 0, 1, 0.94]) 
fig.subplots_adjust(wspace=0) 

plt.show()