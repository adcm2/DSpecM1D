#%%
import numpy as np
import matplotlib.pyplot as plt
import math
import re
from matplotlib import patheffects

path1 = "record_section.out"
d1 = np.loadtxt(path1, delimiter=";")



# initialise plot
maxval = 4
fig, ax = plt.subplots(1,1)
SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20   
maxval = 0 
for i in range(2,88):
    maxval = max(maxval, d1[:,3 * i +2].max())

maxval = maxval/20

# plot
for i in range(0,91):
    maxval = np.max(np.abs(d1[:,3 * i + 2]))/2
    ax.plot(d1[:,3 * i + 2]/maxval + 2 * i, d1[:,0],"k", linewidth=1.6)


# Read in TauP travel time data
data = np.genfromtxt("travel_times.txt", delimiter=";", names=True, dtype=None, encoding='utf-8')
dist_min, dist_max = 0, 180
time_min, time_max = 0, 3500

allphases = ["P", "S", "PcP", "ScS", "SS", "PP", "ScP", "PKP", "SKS", "sP", "sS", "pP", "pS", "pPcP", "sScS", "PKKP", "SKKS","p","s","PcpPcP","ScsScS","sScsScS","PPP","SSS","pPP","sSS","PKiKP","PcpPcpPcP","ScsScsScS","sScsScsScS","Pdiff","Sdiff","PKiKP","SKiKP","sPKiKP","sSKiKP","pPdiff","sSdiff","PS","pPS"]

phases_to_show = ["P", "S", "PcP", "ScS", "p", "s", "sS", "pP", "ScsScS","sScsScS","ScsScsScS","sScsScsScS", "pPcP", "sScS", "Sdiff","SS","pPS","PS","pS"]

phases_to_show = ["s","S", "ScS", "sS", "ScsScS","sScsScS","ScsScsScS","sScsScsScS", "sScS", "Sdiff","SS","sSv410sScS","sSv660sScS","SKS","SSS","SSSS","SSSSS"]
# phases_to_show = ["ScS","sSv410sScS","sSv660sScS"]
# phases_to_show = ["SS"]
# phases_to_show = ["PcS"]
# plt.figure(figsize=(10, 6))
text_style = [patheffects.withStroke(linewidth=3, foreground='white')]

# 3. Filter and plot
for phase_name in phases_to_show:
    # Create a mask for the specific phase
    mask = data['Phase'] == phase_name
    
    # Extract distances and times for this phase
    dist = data['Distance_deg'][mask]
    time = data['TravelTime_s'][mask]
    
    if len(dist) > 0:
        # Since the file is sorted by ray_param, plotting as a line 
        # handles the triplication "bow-tie" correctly.
        ax.plot(dist, time, "r",label=phase_name, linewidth=2.0)

        # FIND THE BEST LABEL POSITION
        # We look for points that are inside the current view limits
        in_view = (dist >= dist_min) & (dist <= dist_max) & \
                  (time >= time_min) & (time <= time_max)
        # idx = np.argmax(dist)
        # plt.text(dist[idx], time[idx], f' {phase_name}', 
        #          color='r', fontsize=10, fontweight='bold',
        #          verticalalignment='center')
        # Get the index of the last point that is still inside the box
        idx = np.where(in_view)[0][-1]
            
        # Add text with path effects (white outline)
        txt = ax.text(dist[idx], time[idx], f" {phase_name}", 
                          color='r', fontsize=11, fontweight='bold',
                          va='center', ha='left')
        txt.set_path_effects(text_style)
    else:
        print(f"Warning: Phase '{phase_name}' not found in file.")
ax.set_ylim(time_min, time_max)
plt.show()