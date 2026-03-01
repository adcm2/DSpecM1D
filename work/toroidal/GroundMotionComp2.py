#%%
import numpy as np
import matplotlib.pyplot as plt
import math

path = "groundresponsel_comp.out"
d = np.loadtxt(path, delimiter=";")
# path1 = "groundresponsel_comp.out"
# d1 = np.loadtxt(path1, delimiter=";")

# initialise plot
maxval = 4
fig, ax = plt.subplots(1,1)
plt.subplots_adjust(wspace=0, hspace=0)
SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20    
# ax[0].plot(d[:,0],d[:,3], "b", linewidth=3)
# ax[0].plot(d[:,0],d[:,6], "r", linewidth=3)
# ax[0].plot(d[:,0],d[:,9], "m--", linewidth=3)
# ax[0].plot(d[:,0],d[:,10], "y--", linewidth=3)

ax.plot(d[:,0],d[:,1], "g", linewidth=3)
ax.plot(d[:,0],d[:,7], "c--", linewidth=3)
ax.plot(d[:,0],d[:,13], "k-.", linewidth=3)
ax.legend(['PREM','Perturbed PREM', 'Unperturbed basis'])

# ax.legend(['Rel error theta'], loc='upper right', fontsize=BIGGER_SIZE)

ax.set_xlabel('Frequency (mHz)', fontsize=BIGGER_SIZE)
ax.set_ylabel('Amplitude', fontsize=BIGGER_SIZE)      
plt.show()