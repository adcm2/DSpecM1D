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

# Filter for frequency range of interest (0.0 to 5.0 mHz)
fmin, fmax = 0.0, 5.0
idx_min = np.searchsorted(dmf[:, 0], fmin)
idx_max = np.searchsorted(dmf[:, 0], fmax, side='right')
print(f"Plotting frequency range: {dmf[idx_min, 0]:.2f} mHz to {dmf[idx_max,0]:.2f} mHz")

time_vector = dmf[idx_min:idx_max, 0]

# --- 4. Extract Data Columns ---
trans_z = dmf[idx_min:idx_max, 3]
yspec_z = dmf[idx_min:idx_max, 12]

# --- 5. Calculate Statistics ---
norm_z = np.max(np.abs(yspec_z)) * 1.01  # Add small buffer

# Calculate relative difference
yspec_diff = np.abs(yspec_z - trans_z) / norm_z * 100
yspec_av_diff = np.mean(yspec_diff)

print(f"Average relative difference (Z): {yspec_av_diff:.4f} %")
print(f"Maximum relative difference (Z): {np.max(yspec_diff):.4f} %")

# --- 6. Plot Data ---
# Add small offset (1e-3) to avoid log scale issues if used later, or just visual base
ax_data.plot(time_vector, yspec_z/norm_z + 1e-3, "b", linewidth=lwidth, label='YSpec')
ax_data.plot(time_vector, trans_z/norm_z + 1e-3, "r--", linewidth=lwidth, label='DSpecM1D')

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

# Plot data on inset
ax_ins.plot(zoom_x, zoom_y_yspec, "b", linewidth=lwidth)
ax_ins.plot(zoom_x, zoom_y_trans, "r--", linewidth=lwidth)

# Calculate max value within the zoom range for the y-limit
zoom_ymax = max(np.max(zoom_y_yspec), np.max(zoom_y_trans)) * 1.07
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

names = [r'${}_0$S${}_2$', r'${}_2$S${}_1$', r'${}_0$S${}_3$']
xvalues = [0.3108, 0.4063, 0.4713]
yvalues = [0.0832, 0.0096, 0.0556]

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

# Plot data on second inset
ax_ins2.plot(zoom2_x, zoom2_y_yspec, "b", linewidth=lwidth)
ax_ins2.plot(zoom2_x, zoom2_y_trans, "r--", linewidth=lwidth)

# Calculate max for y-limit
zoom2_ymax = max(np.max(zoom2_y_yspec), np.max(zoom2_y_trans)) * 1.07
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

names2 = [r'${}_4$S${}_4$', r'${}_1$S${}_{11}$', r'${}_5$S${}_4$',r'${}_2$S${}_{10}$', r'${}_6$S${}_2$', r'${}_0$S${}_{16}$',r'${}_2$S${}_0$',r'${}_7$S${}_2$']
xvalues2 = [2.292460, 2.359480, 2.387861, 2.425544, 2.460181, 2.472586,2.513545 ,2.532304]
yvalues2 = [0.081, 0.57, 0.12, 0.35, 0.005, 0.4, 0.007, 0.026] 

for i, name in enumerate(names2):
    ax_ins2.axvline(xvalues2[i], ymin=0,ymax=(yvalues2[i] + 0.01)/y2_2, color='black', linestyle='--', linewidth=lwidth2)
    if i in [0,1,2,3,5,6,7]:  # For the more prominent peaks, place label above the line
        ax_ins2.text(xvalues2[i], yvalues2[i] + 0.01, name, fontsize=14, fontweight='bold', ha='center')
    else:  # For smaller peaks, place label slightly above the line to avoid overlap
        ax_ins2.text(xvalues2[i]-0.01, yvalues2[i] + 0.015, name, fontsize=14, fontweight='bold', ha='center')


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