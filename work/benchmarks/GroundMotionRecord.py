#%%
import numpy as np
import matplotlib.pyplot as plt
import math
import re

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
for i in range(0,90):
    maxval = np.max(np.abs(d1[:,3 * i + 2]))/2
    ax.plot( d1[:,3 * i + 2]/maxval + 2 * i, d1[:,0]/60,"k", linewidth=1.6)
ax.set_ylim(0, 50)
plt.show()