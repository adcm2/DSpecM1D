#%%
import numpy as np
import matplotlib.pyplot as plt
import math

# path = "fullresponse_single_w.out"
# d = np.loadtxt(path, delimiter=";")

path1 = "full_solution_w.out"
d1 = np.loadtxt(path1, delimiter=";")

# initialise plot
maxval = 4
fig, ax = plt.subplots(1,3)
plt.subplots_adjust(wspace=0, hspace=0)
SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20    
ax[0].plot(d1[:,0],d1[:,1], "b", linewidth=3)
ax[0].plot(d1[:,0],d1[:,10], "r--", linewidth=3)
ax[0].plot(d1[:,0],d1[:,19], "g--", linewidth=3)
ax[0].plot(d1[:,0],d1[:,28], "k--", linewidth=3)
ax[0].legend(['UD','U','VD','V'], loc='upper right', fontsize=BIGGER_SIZE)
ax[0].set_xlim(5250,6371)

ax[1].plot(d1[:,0],d1[:,4], "r", linewidth=3)
ax[1].plot(d1[:,0],d1[:,13], "b", linewidth=3)
ax[1].plot(d1[:,0],d1[:,22], "g", linewidth=3)
ax[1].plot(d1[:,0],d1[:,31], "k", linewidth=3)
ax[1].set_xlim(5250,6371)


ax[2].plot(d1[:,0],d1[:,7], "g", linewidth=3)
ax[2].plot(d1[:,0],d1[:,16], "b", linewidth=3)
ax[2].plot(d1[:,0],d1[:,25], "r", linewidth=3)
ax[2].plot(d1[:,0],d1[:,34], "k", linewidth=3)
ax[2].set_xlim(5250,6371)

# ax[0].plot(d[:,0],d[:,1], "b", linewidth=3)
# ax[0].plot(d[:,0],d[:,2], "r", linewidth=3)
# ax[1].plot(d[:,0],d[:,3], "b", linewidth=3)
# ax[1].plot(d[:,0],d[:,4], "r", linewidth=3)

# ax[2].plot(d[:,0],d[:,5], "b", linewidth=3)
# ax[2].plot(d[:,0],d[:,6], "r", linewidth=3)
# ax.plot(d[:,0],d[:,12], "r", linewidth=3)
plt.show()