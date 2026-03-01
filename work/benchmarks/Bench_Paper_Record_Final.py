#%%
import numpy as np
import matplotlib.pyplot as plt
import math
import re
from matplotlib import patheffects

# --- 1. Load Data ---
path1 = "record_section.out"
try:
    d1 = np.loadtxt(path1, delimiter=";")
except FileNotFoundError:
    print(f"Error: Data file {path1} not found")
    exit()

# Read in TauP travel time data
try:
    data = np.genfromtxt("travel_times.txt", delimiter=";", names=True, dtype=None, encoding='utf-8')
except FileNotFoundError:
    print(f"Error: Data file travel_times.txt not found")
    exit()

# --- 2. Setup Figure ---
fig, ax = plt.subplot_mosaic([['left', 'right']], figsize=(14, 10), constrained_layout=True)
ax_p, ax_s = ax['left'], ax['right']

SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 24

# --- 3. Plot Waveforms ---
# Left Plot (Vertical Component / P-waves area)
for i in range(0, 91):
    maxvalp = np.max(np.abs(d1[:, 3 * i + 1])) / 2
    if maxvalp > 0:
        ax_p.plot(d1[:, 3 * i + 1] / maxvalp + 2 * i, d1[:, 0], "k", linewidth=1.1)

# Right Plot (Horizontal Component / S-waves area)
for i in range(0, 91):
    maxvals = np.max(np.abs(d1[:, 3 * i + 2])) / 2
    if maxvals > 0:
        ax_s.plot(d1[:, 3 * i + 2] / maxvals + 2 * i, d1[:, 0], "k", linewidth=1.1)

# --- 4. Plot Travel Time Curves ---
dist_min, dist_max = 0, 180
time_min, time_max = 0, 3600

phases_s = ["s", "S", "ScS", "sS", "ScsScS", "sScsScS", "ScsScsScS", "sScsScsScS", "sScS", "Sdiff", "SS", "sSv410sScS", "sSv660sScS", "SKS", "SSS", "SSSS", "SSSSS"]
phases_p = ["p", "P", "PcP", "Pdiff", "pP", "PKP", "PP", "pPP","PKIKP","pPPP","PPP","PcPPcP"]

text_style = [patheffects.withStroke(linewidth=3, foreground='white')]

def plot_phases(axes, phase_list, data_struct):
    for phase_name in phase_list:
        mask = data_struct['Phase'] == phase_name
        dist = data_struct['Distance_deg'][mask]
        time = data_struct['TravelTime_s'][mask]
        
        if len(dist) > 0:
            axes.plot(dist, time, "r", label=phase_name, linewidth=2.0)
            
            # Label positioning
            in_view = (dist >= dist_min) & (dist <= dist_max) & \
                      (time >= time_min) & (time <= time_max)
            
            if np.any(in_view):
                idx_candidates = np.where(in_view)[0]
                idx = idx_candidates[-1] # Pick the last point in view
                
                txt = axes.text(dist[idx], time[idx], f" {phase_name}", 
                                color='r', fontsize=14, fontweight='bold',
                                va='center', ha='left')
                txt.set_path_effects(text_style)

# Plot P phases on left, S phases on right
plot_phases(ax_p, phases_p, data)
plot_phases(ax_s, phases_s, data)

# --- 5. Formatting ---
xticks_val = [0, 45, 90, 135, 180]
yticks_val = [0, 600, 1200, 1800, 2400, 3000, 3600]
spine_width = 2.0

for target_ax, title in zip([ax_p, ax_s], ['Vertical (Z)', 'North (T)']):
    target_ax.set_ylim(time_max, time_min) # Reversed Y axis for record sections usually
    target_ax.set_ylim(time_min, time_max) # Standard Y axis direction
    
    target_ax.set_xticks(xticks_val)
    target_ax.set_yticks(yticks_val)
    
    target_ax.set_xlim(dist_min, dist_max)
    target_ax.set_ylim(time_min, time_max)
    
    target_ax.set_xlabel('Distance (deg)', fontsize=BIGGER_SIZE)
    target_ax.tick_params(axis='both', which='major', labelsize=MEDIUM_SIZE, width=spine_width, length=6)
    
    # Thicken spines
    for spine in target_ax.spines.values():
        spine.set_linewidth(spine_width)
        
    target_ax.set_title(title, fontsize=BIGGER_SIZE, fontweight='bold')

ax_p.set_ylabel('Time (s)', fontsize=BIGGER_SIZE)

# Ensure y-labels only on left plot if desired, or keep both
# ax_s.set_yticklabels([]) # Uncomment to hide y-labels on second plot

plt.show()