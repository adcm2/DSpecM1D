#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch

# --- 1. Load Data ---
try:
    path_mf = "../outputs/ex1_w.out"
    dmf = np.loadtxt(path_mf, delimiter=";")
except FileNotFoundError:
    print(f"Error: Data file not found at '{path_mf}'")
    exit()

# --- 2. Setup Figure ---
fig, ax_data = plt.subplots(1, 1, figsize=(14, 12), sharex=True)
plt.style.use('seaborn-v0_8-whitegrid')

# --- 3. Define Plotting Parameters ---
lwidth = 1.3
lwidth2 = 2.0
M_SIZE = 15
BIGGER_SIZE = 20
offset = 10  # The amount of offset in points

# Colorblind-friendly palette (Okabe-Ito)
COLORS = {
    'yspec': '#0072B2',   # blue
    'specnm': '#009E73',  # bluish green
    'dspecm': '#D55E00',  # vermillion
    'text': '#111111'
}

# Filter for frequency range of interest (0.0 to 5.0 mHz)
fmin, fmax = 0.0, 5.0
idx_min = np.searchsorted(dmf[:, 0], fmin)
idx_max = np.searchsorted(dmf[:, 0], fmax, side='right')
print(f"Plotting frequency range: {dmf[idx_min, 0]:.2f} mHz to {dmf[idx_max,0]:.2f} mHz")

time_vector = dmf[idx_min:idx_max, 0]

# --- 4. Extract Data Columns ---
trans_z = dmf[idx_min:idx_max, 3]
yspec_z = dmf[idx_min:idx_max, 12]
# If a third comparison block exists in ex1_w.out, use its Z-abs column.
# Otherwise, fall back to the second block to keep plotting compatible.
specnm_col = 21 
# if dmf.shape[1] > 21 else 12
specnm_z = dmf[idx_min:idx_max, specnm_col]

# --- 5. Calculate Statistics ---
norm_z = np.max(np.abs(yspec_z)) * 1.01  # Add small buffer

# Calculate relative difference
yspec_diff = np.abs(yspec_z - trans_z) / norm_z * 100
yspec_av_diff = np.mean(yspec_diff)
specnm_diff = np.abs(trans_z - specnm_z) / norm_z * 100
specnm_av_diff = np.mean(specnm_diff)

print(f"Average relative difference with yspec  (Z): {yspec_av_diff:.4f} %")
print(f"Maximum relative difference with yspec  (Z): {np.max(yspec_diff):.4f} %")
print(f"Average relative difference with specnm (Z): {specnm_av_diff:.4f} %")
print(f"Maximum relative difference with specnm (Z): {np.max(specnm_diff):.4f} %")

# --- 6. Plot Data ---
# Add small offset (1e-3) to avoid log scale issues if used later, or just visual base
ax_data.plot(time_vector, yspec_z/norm_z + 1e-3, color=COLORS['yspec'], linestyle='-', linewidth=lwidth, label='YSpec')
ax_data.plot(time_vector, trans_z/norm_z + 1e-3, color=COLORS['dspecm'], linestyle='--', linewidth=lwidth, label='DSpecM1D')
ax_data.plot(time_vector, specnm_z/norm_z + 1e-3, color=COLORS['specnm'], linestyle='-.', linewidth=lwidth, label='SpecNM')


def get_modes_in_range(fmin, fmax):
    mode_files = [
        "/home/adcm2/space/c++/mineos/DEMO/MYEX/lf_prem_R",
        "/home/adcm2/space/c++/mineos/DEMO/MYEX/lf_prem_S"
    ]
    modes = []
    for filepath in mode_files:
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
                start_idx = 0
                for i, line in enumerate(lines):
                    if "mode" in line and "phs vel" in line and "w(mhz)" in line:
                        start_idx = i + 2
                        break
                if start_idx > 0:
                    for line in lines[start_idx:]:
                        if not line.strip():
                            continue
                        parts = line.split()
                        if len(parts) >= 5:
                            try:
                                n = parts[0]
                                type_mode = parts[1].upper()
                                l = parts[2]
                                freq = float(parts[4])
                                if fmin <= freq <= fmax:
                                    name = f'${{_{n}}}${type_mode}${{_{{{l}}}}}$'
                                    modes.append((name, freq))
                            except ValueError:
                                pass
        except FileNotFoundError:
            pass
    modes.sort(key=lambda x: x[1])
    return [m[0] for m in modes], [m[1] for m in modes]

# --- 7. Inset Plot (Zoomed Region) ---
# Create inset axes in the TOP LEFT corner
# [left, bottom, width, height] in normalized (0, 1) units of the main axes
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
# ax_ins = ax_data.inset_axes([0.02, 0.4, 0.35, 0.5]) 
ax_ins = fig.add_axes([0.07, 0.45, 0.15, 0.45])  # Alternative way to create inset axes
# ax_ins = zoomed_inset_axes(ax_data, zoom=2.5, loc=1)  # loc=1 is upper right, but we will adjust it to top left
# Define NEW zoom range (0.2 -> 0.5 mHz)
zoom_fmin, zoom_fmax = 0.28, 0.5
x1, x2 = zoom_fmin, zoom_fmax

idx_zmin = np.searchsorted(time_vector, zoom_fmin)
idx_zmax = np.searchsorted(time_vector, zoom_fmax, side='right')

# Get the data relative to the zoom window for plotting and scaling
zoom_x = time_vector[idx_zmin:idx_zmax]
zoom_y_yspec = yspec_z[idx_zmin:idx_zmax]/norm_z + 1e-3
zoom_y_trans = trans_z[idx_zmin:idx_zmax]/norm_z + 1e-3
zoom_y_specnm = specnm_z[idx_zmin:idx_zmax]/norm_z + 1e-3

# Plot data on inset
ax_ins.plot(zoom_x, zoom_y_yspec, color=COLORS['yspec'], linestyle='-', linewidth=lwidth)
ax_ins.plot(zoom_x, zoom_y_trans, color=COLORS['dspecm'], linestyle='--', linewidth=lwidth)
ax_ins.plot(zoom_x, zoom_y_specnm, color=COLORS['specnm'], linestyle='-.', linewidth=lwidth)

# Calculate max value within the zoom range for the y-limit
zoom_ymax = max(np.max(zoom_y_yspec), np.max(zoom_y_trans), np.max(zoom_y_specnm)) * 1.07
y1,y2 = 0, zoom_ymax  

# Style Inset
ax_ins.set_xlim(zoom_fmin, zoom_fmax)
ax_ins.set_ylim(0, zoom_ymax) 
ax_ins.tick_params(axis='x', which='major', labelsize=12, labelcolor='black', colors='black')
ax_ins.set_yticks([]) # Switch off tick values on y axis

ax_ins.grid(False) # Turn off grid for the inset

# Make inset spines (borders) black
for spine in ax_ins.spines.values():
    spine.set_edgecolor('black')
    spine.set_linewidth(lwidth2)


con = ConnectionPatch(
    xyA=(x1, y2), coordsA=ax_data.transData,   # Start: Top-left of zoom box
    xyB=(0, 0),   coordsB=ax_ins.transAxes, # End: Bottom-left of inset
    color="black", linestyle="-",linewidth=lwidth2
)
fig.add_artist(con)
con1 = ConnectionPatch(
    xyA=(x2, y2), coordsA=ax_data.transData,   # Start: Top-left of zoom box
    xyB=(1, 0),   coordsB=ax_ins.transAxes, # End: Bottom-left of inset
    color="black", linestyle="-",linewidth=lwidth2
)
fig.add_artist(con1)

rect = plt.Rectangle((x1, y1), x2-x1, y2-y1, facecolor='none', edgecolor='black', linewidth=lwidth2)
ax_data.add_patch(rect)


names, xvalues = get_modes_in_range(zoom_fmin, zoom_fmax)
print(names)
print(xvalues)
yvalues = []
for x in xvalues:
    idx = np.searchsorted(zoom_x, x)
    if idx < len(zoom_x):
        window = 3
        start_w = max(0, idx - window)
        end_w = min(len(zoom_x), idx + window + 1)
        y_val1 = np.max(zoom_y_yspec[start_w:end_w]) if len(zoom_y_yspec[start_w:end_w]) > 0 else 0
        y_val2 = np.max(zoom_y_trans[start_w:end_w]) if len(zoom_y_trans[start_w:end_w]) > 0 else 0
        y_val3 = np.max(zoom_y_specnm[start_w:end_w]) if len(zoom_y_specnm[start_w:end_w]) > 0 else 0
        yvalues.append(max(y_val1, y_val2, y_val3))
    else:
        yvalues.append(0.05)


for i, name in enumerate(names):
    ax_ins.axvline(xvalues[i], ymin=0,ymax=(yvalues[i]+0.002)/y2, color='black', linestyle='--', linewidth=lwidth2)
    ax_ins.text(xvalues[i], yvalues[i] + 0.002, name, fontsize=14, fontweight='bold', ha='center')


# --- 7b. Second Inset Plot (2.25 - 2.55 mHz) ---
ax_ins2 = fig.add_axes([0.25, 0.45, 0.15, 0.45])  # Position to the right of the first inset

zoom2_fmin, zoom2_fmax = 2.27, 2.55
x1_2, x2_2 = zoom2_fmin, zoom2_fmax

idx_zmin2 = np.searchsorted(time_vector, zoom2_fmin)
idx_zmax2 = np.searchsorted(time_vector, zoom2_fmax, side='right')

# Get data for second zoom
zoom2_x = time_vector[idx_zmin2:idx_zmax2]
zoom2_y_yspec = yspec_z[idx_zmin2:idx_zmax2]/norm_z + 1e-3
zoom2_y_trans = trans_z[idx_zmin2:idx_zmax2]/norm_z + 1e-3
zoom2_y_specnm = specnm_z[idx_zmin2:idx_zmax2]/norm_z + 1e-3

# Plot data on second inset
ax_ins2.plot(zoom2_x, zoom2_y_yspec, color=COLORS['yspec'], linestyle='-', linewidth=lwidth)
ax_ins2.plot(zoom2_x, zoom2_y_trans, color=COLORS['dspecm'], linestyle='--', linewidth=lwidth)
ax_ins2.plot(zoom2_x, zoom2_y_specnm, color=COLORS['specnm'], linestyle='-.', linewidth=lwidth)

# Calculate max for y-limit
zoom2_ymax = max(np.max(zoom2_y_yspec), np.max(zoom2_y_trans), np.max(zoom2_y_specnm)) * 1.07
y1_2, y2_2 = 0, zoom2_ymax

# Style Second Inset
ax_ins2.set_xlim(zoom2_fmin, zoom2_fmax)
ax_ins2.set_ylim(0, zoom2_ymax)
ax_ins2.tick_params(axis='x', which='major', labelsize=12, labelcolor='black', colors='black')
ax_ins2.set_yticks([]) 
ax_ins2.grid(False)

for spine in ax_ins2.spines.values():
    spine.set_edgecolor('black')
    spine.set_linewidth(lwidth2)

# Connection lines for second inset
con2_a = ConnectionPatch(
    xyA=(x1_2, y2_2), coordsA=ax_data.transData, 
    xyB=(1, 0), coordsB=ax_ins2.transAxes, 
    color="black", linestyle="-", linewidth=lwidth2
)
fig.add_artist(con2_a)

con2_b = ConnectionPatch(
    xyA=(x2_2, y2_2), coordsA=ax_data.transData, 
    xyB=(1, 1), coordsB=ax_ins2.transAxes, 
    color="black", linestyle="-", linewidth=lwidth2
)
fig.add_artist(con2_b)

# Rectangle on main plot for second range
rect2 = plt.Rectangle((x1_2, y1_2), x2_2-x1_2, y2_2-y1_2, facecolor='none', edgecolor='black', linewidth=lwidth2)
ax_data.add_patch(rect2)


names2, xvalues2 = get_modes_in_range(zoom2_fmin, zoom2_fmax)
print(names2)
print(xvalues2)
yvalues2 = []
for x in xvalues2:
    idx = np.searchsorted(zoom2_x, x)
    if idx < len(zoom2_x):
        window = 3
        start_w = max(0, idx - window)
        end_w = min(len(zoom2_x), idx + window + 1)
        y_val1 = np.max(zoom2_y_yspec[start_w:end_w]) if len(zoom2_y_yspec[start_w:end_w]) > 0 else 0
        y_val2 = np.max(zoom2_y_trans[start_w:end_w]) if len(zoom2_y_trans[start_w:end_w]) > 0 else 0
        y_val3 = np.max(zoom2_y_specnm[start_w:end_w]) if len(zoom2_y_specnm[start_w:end_w]) > 0 else 0
        yvalues2.append(max(y_val1, y_val2, y_val3))
    else:
        yvalues2.append(0.05)
 

for i, name in enumerate(names2):
    ax_ins2.axvline(xvalues2[i], ymin=0,ymax=(yvalues2[i] + 0.01)/y2_2, color='black', linestyle='--', linewidth=lwidth2)

    if i % 2 == 0:  # Alternate label placement height for clear reading
        ax_ins2.text(xvalues2[i], yvalues2[i] + 0.01, name, fontsize=14, fontweight='bold', ha='center')
    else:  
        ax_ins2.text(xvalues2[i], yvalues2[i] + 0.02, name, fontsize=14, fontweight='bold', ha='center')


# --- 8. Styling ---
# Axis Labels & Title
ax_data.set_xlabel('Frequency (mHz)', fontsize=BIGGER_SIZE)
ax_data.set_ylabel('Z displacement (N.D.)', fontsize=BIGGER_SIZE)

# Legend
ax_data.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=3, 
               fancybox=True, shadow=True, fontsize=BIGGER_SIZE, 
               frameon=True, edgecolor='black')

# Spines (Borders)
ax_data.spines['top'].set_visible(False)
ax_data.spines['right'].set_visible(False)
ax_data.spines['left'].set_position(('outward', offset))
ax_data.spines['left'].set_linewidth(1.5)
ax_data.spines['bottom'].set_position(('outward', offset))
ax_data.spines['bottom'].set_linewidth(1.5)

# Ticks
ax_data.tick_params(axis='both', which='major', labelsize=M_SIZE, 
                    length=10, width=1.5, labelcolor='black')

# Limits
ax_data.set_xlim([fmin, fmax])
z_axlim = 1.0
ax_data.set_ylim(0, z_axlim)
ax_data.set_yticks([0, z_axlim])

plt.tight_layout()
plt.show()