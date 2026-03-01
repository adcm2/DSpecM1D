#%%
import numpy as np
import matplotlib.pyplot as plt
import math
import re




# initialise plot
maxval = 4
fig, ax = plt.subplots(1,2)
SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20   

for i in range(1,60):
    path = "rad_start_" + str(i+1) + ".out"
    d = np.loadtxt(path, delimiter=";")
    ax[0].plot(d[:,0], d[:,1], linewidth=2)
    ax[0].legend(['l=2', 'l=3'], fontsize=14)
    dd = np.diff(d[:,1])/np.diff(d[:,0])
    # for j in range(len(dd)):
    #     if ()
    #     dd[j] = dd[j]/(d[j+1,0]-d[j,0])

    # dm = np.max(np.abs(dd))
    ax[1].plot(d[1:,0], dd, linewidth=2)


# ax.legend(['l=2', 'l=3', 'l=4', 'l=5', 'l=6', 'l=7', 'l=8', 'l=9', 'l=10'], fontsize=14)

plt.show()