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
    path_mf = "tidal_radial_response0.010000.out"
    dmf = np.loadtxt(path_mf, delimiter=";")
    path_mf2 = "tidal_radial_response0.001000.out"
    dmf2 = np.loadtxt(path_mf2, delimiter=";")
    path_mf3 = "tidal_radial_response0.000100.out"
    dmf3 = np.loadtxt(path_mf3, delimiter=";")

    path_freqs = "tidal_frequency_0.010000.out"
    dfreqs = np.loadtxt(path_freqs, delimiter=";")
    # Assuming freqs are same structure for all, or load distict if needed
    
    path_n2 = "tidal_N2_0.010000.out"
    dn2 = np.loadtxt(path_n2, delimiter=";")
    # ... assuming N2 is same or similar for display purposes ...
    
except FileNotFoundError:
    print(f"Error: Data file not found")
    exit()

# --- 2. Setup a Cleaner Figure ---
nfreqs = dfreqs.shape[0]  # Number of frequency columns
if(nfreqs > 6):
    nfreqs = 6  # Limit to first 6 frequencies for plotting

# Create a 3-row grid (one for each dataset)
fig, axes = plt.subplots(3, nfreqs + 1, figsize=(14, 15), sharey=True)
fig.subplots_adjust(wspace=0, hspace=0.1) 
plt.style.use('seaborn-v0_8-whitegrid')

# --- 3. Define Plotting Parameters ---
lw = 1.2
M_SIZE = 12
BIGGER_SIZE = 18
time_vector = dmf[:, 0]

# --- Data Column Indices ---
# Give meaningful names to the columns being plotted.
# This makes the code MUCH easier to read and debug.
rad_val = dmf[:, 0] # They likely share the same radial vector if grid is same
rad_val2 = dmf2[:, 0]
rad_val3 = dmf3[:, 0]

low_rad = dn2[0,0]
up_rad = dn2[-1,0]
# List of Roman numerals for labeling
roman_numerals = ['i)', 'ii)', 'iii)', 'iv)', 'v)', 'vi)']
viscosity_labels = ['visc=0.01', 'visc=0.001', 'visc=0.0001']

# Helper function to plot a row
def plot_row(row_idx, data_mf, rad_v, label_prefix=""):
    for i in range(nfreqs):
        ax_data = axes[row_idx, i]
        
        # Calculate max for scaling
        maxvalt = 1.1*max(
            np.max(np.abs(data_mf[:, 6 * i + 1])), 
            np.max(np.abs(data_mf[:, 6 * i + 4]))
        )
        if maxvalt == 0: maxvalt = 1.0 # Avoid singular limits
        
        ax_data.plot(data_mf[:, 6 * i + 1], rad_v, "b", linewidth=lw, label='U')
        ax_data.plot(data_mf[:, 6 * i + 4], rad_v, "r", linewidth=lw, label='V')

        # Labeling
        # label_text = f"{roman_numerals[i]} {dfreqs[i,0]:.1f}"
        label_text = f"{roman_numerals[i]}"
        ax_data.text(0.05, 0.95, label_text, transform=ax_data.transAxes, 
                     fontsize=M_SIZE, verticalalignment='top', fontweight='bold')
        
        # Viscosity label on first column of each row
        # if i == 0:
            # ax_data.text(0.05, 0.85, viscosity_labels[row_idx], transform=ax_data.transAxes,
                        # fontsize=M_SIZE-2, verticalalignment='top', color='black')

        # Lines and Styling
        ax_data.axhline(low_rad, color='k', linestyle='--', linewidth=1)
        ax_data.axhline(up_rad, color='k', linestyle='--', linewidth=1)
        ax_data.axvline(0, color='k', linestyle='--', linewidth=1)
        
        ax_data.set_ylim([0, rad_v[-1]])
        ax_data.set_xlim([-maxvalt, maxvalt])
        
        ax_data.set_xticks([0]) 
        ax_data.set_xticklabels([0])
        
        ax_data.spines['left'].set_linewidth(1.5)
        ax_data.spines['right'].set_linewidth(1.5)
        ax_data.spines['bottom'].set_linewidth(1.5)
        ax_data.spines['top'].set_linewidth(1.5)
        ax_data.tick_params(axis='x', which='major', length=5, width=1.5, labelcolor='black', labelsize=M_SIZE)

        if i > 0:
            ax_data.set_yticklabels([])
            ax_data.set_yticks([])
        if row_idx < 2:
            ax_data.set_xticklabels([])
            ax_data.set_xticks([])

# Plot Row 1 (0.01)
plot_row(0, dmf, rad_val)

# Plot Row 2 (0.001)
plot_row(1, dmf2, rad_val2)

# Plot Row 3 (0.0001)
plot_row(2, dmf3, rad_val3)


# Legend on top row, first plot
ax_legend = axes[0, 0]
legend_properties = {'weight': 'bold', 'size': BIGGER_SIZE}
ax_legend.legend(loc='upper center', bbox_to_anchor=(2.5, 1.25), 
               ncol=4, fancybox=True, shadow=True, 
               frameon=True, edgecolor='black', prop=legend_properties)

# N^2 Profile (on all 3 rows? Or just one? Usually just one is enough, but grid requires filling)
# Let's plot N2 on all 3 rows for consistency or leave blank. Plotting on all for symmetry.
for r in range(3):
    ax_n2 = axes[r, nfreqs]
    ax_n2.plot(dn2[:,1], dn2[:,0],"g", linewidth=lw)
    ax_n2.axvline(0, color='k', linestyle='--', linewidth=1)
    ax_n2.axhline(low_rad, color='k', linestyle='--', linewidth=1)
    ax_n2.axhline(up_rad, color='k', linestyle='--', linewidth=1)
    ax_n2.set_xticks([0])
    ax_n2.set_xticklabels([0])
    ax_n2.set_yticklabels([])
    ax_n2.set_yticks([])
    
    # Styling
    ax_n2.spines['left'].set_linewidth(1.5)
    ax_n2.spines['right'].set_linewidth(1.5)
    ax_n2.spines['bottom'].set_linewidth(1.5)
    ax_n2.spines['top'].set_linewidth(1.5)
    
    ax_n2.text(0.05, 0.95, "v)", transform=ax_n2.transAxes, 
                 fontsize=M_SIZE, verticalalignment='top', fontweight='bold')
    ax_n2.tick_params(axis='x', which='major', length=5, width=1.2,labelcolor='black', labelsize=M_SIZE)
    if r < 2:
        ax_n2.set_xticklabels([])
        ax_n2.set_xticks([])

# Use tight_layout first to crop outer whitespace
# Reduced pad from 0.5 to 0.1 to minimize outer border
plt.tight_layout(pad=0.1, w_pad=0.0, h_pad=0.0)

# Then adjust internal spacing
fig.subplots_adjust(wspace=0, hspace=0.05) 
plt.show()