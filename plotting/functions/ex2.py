#%%
import numpy as np
import matplotlib.pyplot as plt
import math
import re

# =============================================================================
# 1. HELPER FUNCTIONS & DATA LOADING
# =============================================================================

# This helper function reads the data file.
_float_re = re.compile(r'[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][+-]?\d+)?')
def read_spectra_to_array(path, pad=np.nan):
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

# Load the primary dataset
try:
    path_mf = "../outputs/ex2_t.out"
    dmf = np.loadtxt(path_mf, delimiter=";")
except FileNotFoundError:
    print(f"Error: Data file not found at '{path_mf}'")
    exit()

# =============================================================================
# 2. DATA PROCESSING & PARAMETER DEFINITION
# =============================================================================

# --- Define Plotting Parameters ---
lwidth = 2.0
M_SIZE = 15
BIGGER_SIZE = 20

# --- Extract Data Columns ---
lowidx = 0
upidx = 18001
time_vector = dmf[lowidx:upidx, 0]
dspecm_z, dspecm_n, dspecm_e = dmf[lowidx:upidx, 1], dmf[lowidx:upidx, 2], dmf[lowidx:upidx, 3]
yspec_z, yspec_n, yspec_e = dmf[lowidx:upidx, 4], dmf[lowidx:upidx, 5], dmf[lowidx:upidx, 6]
mineos_z, mineos_n, mineos_e = dmf[lowidx:upidx, 7], dmf[lowidx:upidx, 8], dmf[lowidx:upidx, 9]
specnm_z, specnm_n, specnm_e = dmf[lowidx:upidx, 10], dmf[lowidx:upidx, 11], dmf[lowidx:upidx, 12]

# --- Define Time Limits ---
tlow = time_vector[0]
tup = 18000

# --- Calculate Normalization Factors ---
norm_z = max(np.max(np.abs(yspec_z)),np.max(np.abs(dspecm_z)), np.max(np.abs(mineos_z))) + 1e-12
norm_n = max(np.max(np.abs(yspec_n)),np.max(np.abs(dspecm_n)), np.max(np.abs(mineos_n))) + 1e-12
norm_e = max(np.max(np.abs(yspec_e)),np.max(np.abs(dspecm_e)), np.max(np.abs(mineos_e))) + 1e-12

# --- Calculate Average Differences ---
yspec_diff_z = np.abs(yspec_z - dspecm_z) / norm_z * 100
yspec_av_diff_z = np.mean(yspec_diff_z)
mineos_diff_z = np.abs(mineos_z - dspecm_z) / norm_z * 100
mineos_av_diff_z = np.mean(mineos_diff_z)
specnm_diff_z = np.abs(specnm_z - dspecm_z) / norm_z * 100
specnm_av_diff_z = np.mean(specnm_diff_z)
specnm_yspec_diff_z = np.abs(specnm_z - yspec_z) / norm_z * 100
specnm_yspec_av_diff_z = np.mean(specnm_yspec_diff_z)

yspec_diff_n = np.abs(yspec_n - dspecm_n) / norm_n * 100
yspec_av_diff_n = np.mean(yspec_diff_n)
mineos_diff_n = np.abs(mineos_n - dspecm_n) / norm_n * 100
mineos_av_diff_n = np.mean(mineos_diff_n)
specnm_diff_n = np.abs(specnm_n - dspecm_n) / norm_n * 100
specnm_av_diff_n = np.mean(specnm_diff_n)
specnm_yspec_diff_n = np.abs(specnm_n - yspec_n) / norm_n * 100
specnm_yspec_av_diff_n = np.mean(specnm_yspec_diff_n)

yspec_diff_e = np.abs(yspec_e - dspecm_e) / norm_e * 100
yspec_av_diff_e = np.mean(yspec_diff_e)
mineos_diff_e = np.abs(mineos_e - dspecm_e) / norm_e * 100
mineos_av_diff_e = np.mean(mineos_diff_e)
specnm_diff_e = np.abs(specnm_e - dspecm_e) / norm_e * 100
specnm_av_diff_e = np.mean(specnm_diff_e)
specnm_yspec_diff_e = np.abs(specnm_e - yspec_e) / norm_e * 100
specnm_yspec_av_diff_e = np.mean(specnm_yspec_diff_e)

# print max relative differences for debugging
print(f"Max relative difference for Z: YSpec={np.max(yspec_diff_z):.2f} %, MINEOS={np.max(mineos_diff_z):.2f} %, specnm={np.max(specnm_diff_z):.2f} %")
print(f"Max relative difference for N: YSpec={np.max(yspec_diff_n):.2f} %, MINEOS={np.max(mineos_diff_n):.2f} %, specnm={np.max(specnm_diff_n):.2f} %")
print(f"Max relative difference for E: YSpec={np.max(yspec_diff_e):.2f} %, MINEOS={np.max(mineos_diff_e):.2f} %, specnm={np.max(specnm_diff_e):.2f} %")
print(f"Max relative difference between specnm and yspec for Z: {np.max(specnm_yspec_diff_z):.2f} %")
print(f"Max relative difference between specnm and yspec for N: {np.max(specnm_yspec_diff_n):.2f} %")
print(f"Max relative difference between specnm and yspec for E: {np.max(specnm_yspec_diff_e):.2f} %")
print(f"Average relative difference between specnm and yspec for Z: {specnm_yspec_av_diff_z:.2f} %")
print(f"Average relative difference between specnm and yspec for N: {specnm_yspec_av_diff_n:.2f} %")
print(f"Average relative difference between specnm and yspec for E: {specnm_yspec_av_diff_e:.2f} %")

# scale the data by the normalization factors for better visualization
yspec_z /= norm_z
yspec_n /= norm_n
yspec_e /= norm_e
dspecm_z /= norm_z
dspecm_n /= norm_n
dspecm_e /= norm_e
mineos_z /= norm_z
mineos_n /= norm_n
mineos_e /= norm_e
specnm_z /= norm_z
specnm_n /= norm_n
specnm_e /= norm_e
# =============================================================================
# 3. FIGURE SETUP & PLOTTING
# =============================================================================

# --- Create Figure and Axes ---
fig, axes = plt.subplots(3, 1, figsize=(14, 12), sharex=True)
plt.style.use('seaborn-v0_8-whitegrid')
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)

# --- Plot Z Component (Top) ---
ax_data = axes[0]
ax_data.plot(time_vector, yspec_z/1.01, "b", linewidth=lwidth, label='YSpec')
# ax_data.plot(time_vector, mineos_z/1.01, "g--", linewidth=lwidth, label='MINEOS')
ax_data.plot(time_vector, dspecm_z/1.01, "r-.", linewidth=lwidth, label='DSpecM1D')
ax_data.plot(time_vector, specnm_z/1.01, "m:", linewidth=lwidth, label='SpecNM')


# --- Plot North Component (Middle) ---
ax_data = axes[1]
ax_data.plot(time_vector, yspec_n/1.01, "b", linewidth=lwidth)
# ax_data.plot(time_vector, mineos_n/1.01, "g--", linewidth=lwidth)
ax_data.plot(time_vector, dspecm_n/1.01, "r-.", linewidth=lwidth)
ax_data.plot(time_vector, specnm_n/1.01, "m:", linewidth=lwidth)


# --- Plot East Component (Bottom) ---
ax_data = axes[2]
ax_data.plot(time_vector, yspec_e/1.01, "b", linewidth=lwidth)
# ax_data.plot(time_vector, mineos_e/1.01, "g--", linewidth=lwidth)
ax_data.plot(time_vector, dspecm_e/1.01, "r-.", linewidth=lwidth)
ax_data.plot(time_vector, specnm_e/1.01, "m:", linewidth=lwidth)


# =============================================================================
# 4. AESTHETIC ADJUSTMENTS & STYLING
# =============================================================================

# --- Offset Bottom Axes ---
offset = 10  # The amount of offset in points

# --- Style Top Plot (Z-Component) ---
ax_data = axes[0]
ax_data.legend(loc='upper center', bbox_to_anchor=(0.5, 1.3), ncol=3, fancybox=True, shadow=True, fontsize=BIGGER_SIZE, frameon=True, edgecolor='black')
ax_data.spines['top'].set_visible(False)
ax_data.spines['bottom'].set_visible(False)
ax_data.spines['right'].set_visible(False)
ax_data.spines['left'].set_position(('outward', offset))
ax_data.spines['left'].set_linewidth(1.5)
ax_data.tick_params(axis='x', which='both', bottom=False, top=False)
ax_data.tick_params(axis='y', which='major', length=10, width=1.5)
# ax_data.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax_data.tick_params(axis='both', which='major', labelcolor='black',labelsize=M_SIZE)
z_mval = max(np.max(np.abs(yspec_z)),np.max(np.abs(dspecm_z)), np.max(np.abs(mineos_z)))
z_axlim = 1.0 * z_mval
ax_data.set_ylim(-z_axlim, z_axlim)
extraticks = [-z_axlim, 0,  z_axlim]
ax_data.set_yticks(extraticks)
# ax_data.offset_text.set_visible(False)  # Hide the default offset text

# --- Style Middle Plot (N-Component) ---
ax_data = axes[1]
ax_data.spines['top'].set_visible(False)
ax_data.spines['bottom'].set_visible(False)
ax_data.spines['right'].set_visible(False)
ax_data.spines['left'].set_position(('outward', offset))
ax_data.spines['left'].set_linewidth(1.5)
ax_data.tick_params(axis='x', which='both', bottom=False, top=False)
ax_data.tick_params(axis='y', which='major', length=10, width=1.5)
# ax_data.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax_data.tick_params(axis='both', which='major', labelcolor='black',labelsize=M_SIZE)
n_mval = max(np.max(np.abs(yspec_n)),np.max(np.abs(dspecm_n)), np.max(np.abs(mineos_n)))
n_axlim = 1.0 * n_mval
ax_data.set_ylim(-n_axlim, n_axlim)
extraticks = [-n_axlim,  0,  n_axlim]
ax_data.set_yticks(extraticks)

# --- Style Bottom Plot (E-Component) ---

ax_data = axes[2]
ax_data.spines['top'].set_visible(False)
ax_data.spines['right'].set_visible(False)
ax_data.spines['bottom'].set_position(('outward', offset))
ax_data.spines['bottom'].set_linewidth(1.5)
ax_data.spines['left'].set_position(('outward', offset))
ax_data.spines['left'].set_linewidth(1.5)
ax_data.tick_params(axis='both', which='major', length=10, width=1.5)
# ax_data.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))
ax_data.tick_params(axis='both', which='major', labelcolor='black',labelsize=M_SIZE)
e_mval = max(np.max(np.abs(yspec_e)),np.max(np.abs(dspecm_e)), np.max(np.abs(mineos_e)))
e_axlim = 1.0 * e_mval
ax_data.set_ylim(-e_axlim, e_axlim)
extraticks = [-e_axlim, 0,  e_axlim]
ax_data.set_yticks(extraticks)
xticks = [0, 3000, 6000, 9000, 12000, 15000, 18000]
ax_data.set_xticks(xticks)
ax_data.get_xticklabels()[-1].set_horizontalalignment('right')  # Align the last xtick label to the right
ax_data.get_xticklabels()[-1].set_visible(False)

# --- Set Shared X-Axis Properties ---
axes[2].set_xlabel('Time (s)', fontsize=BIGGER_SIZE)
for ax in axes:
    ax.set_xlim(tlow, tup)

# --- Add Text Annotations ---
txval = tlow 
axes[0].text(txval, 0.8 * z_mval, f"{yspec_av_diff_z:.2f} %", fontsize=M_SIZE, color='black')
axes[0].text(txval, 0.6 * z_mval, f"{specnm_av_diff_z:.2f} %", fontsize=M_SIZE, color='green')
axes[0].text(tup, -z_mval, "Z displacement", fontsize=BIGGER_SIZE, color='black', ha='right')

axes[1].text(txval, 0.8 * n_mval, f"{yspec_av_diff_n:.2f} %", fontsize=M_SIZE, color='black')
axes[1].text(txval, 0.6 * n_mval, f"{specnm_av_diff_n:.2f} %", fontsize=M_SIZE, color='green')
axes[1].text(tup, -n_mval, "N displacement", fontsize=BIGGER_SIZE, color='black', ha='right')

axes[2].text(txval, 0.8 * e_mval, f"{yspec_av_diff_e:.2f} %", fontsize=M_SIZE, color='black')
axes[2].text(txval, 0.6 * e_mval, f"{specnm_av_diff_e:.2f} %", fontsize=M_SIZE, color='green')
axes[2].text(tup, - e_mval, "E displacement", fontsize=BIGGER_SIZE, color='black', ha='right')



# =============================================================================
# 5. FINAL LAYOUT & DISPLAY
# =============================================================================
plt.tight_layout(rect=[0, 0, 1, 0.94]) # Adjust rect to make space for suptitle
plt.show()
