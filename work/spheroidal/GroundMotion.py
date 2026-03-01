#%%
import numpy as np
import matplotlib.pyplot as plt
import math

path = "groundresponse.out"
d = np.loadtxt(path, delimiter=";")
# path1 = "groundresponsel_filt.out"
# d1 = np.loadtxt(path1, delimiter=";")

# initialise plot
maxval = 4
fig, ax = plt.subplots(1)
plt.subplots_adjust(wspace=0, hspace=0)
SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20    

# ax[0].plot(d[:,0],d[:,1], "b", linewidth=3)
# ax[0].plot(d[:,0],d[:,2], "r", linewidth=3)
ax.plot(d[:,0],d[:,3],"k", linewidth=3)

# ax[1].plot(d[:,0],d[:,4], "b", linewidth=3)
# ax[1].plot(d[:,0],d[:,5], "r", linewidth=3)
# ax[1].plot(d[:,0],d[:,6],"k", linewidth=3)

# ax[2].plot(d[:,0],d[:,7], "b", linewidth=3)
# ax[2].plot(d[:,0],d[:,8], "r", linewidth=3)
# ax[2].plot(d[:,0],d[:,9],"k", linewidth=3)

ax.set_title('Z amplitude', fontsize=BIGGER_SIZE)
# ax[1].set_title('Theta amplitude', fontsize=BIGGER_SIZE)
# ax[2].set_title('Phi amplitude', fontsize=BIGGER_SIZE)

ax.set_xlabel('Frequency (mHz)', fontsize=BIGGER_SIZE)
# ax[1].set_xlabel('Frequency (mHz)', fontsize=BIGGER_SIZE)
# ax[2].set_xlabel('Frequency (mHz)', fontsize=BIGGER_SIZE)

# ax[0].legend(['Real Z', 'Imag Z', 'Abs Z'], loc='upper right', fontsize=BIGGER_SIZE)
# ax[1].legend(['Real Theta', 'Imag Theta', 'Abs Theta'], loc='upper right', fontsize=BIGGER_SIZE)
# ax[2].legend(['Real Phi', 'Imag Phi', 'Abs Phi'], loc='upper right', fontsize=BIGGER_SIZE)

# ax[0].set_xlabel('Frequency (mHz)', fontsize=BIGGER_SIZE)
# ax[0].set_ylabel('Amplitude', fontsize=BIGGER_SIZE)      
plt.show()