#%%
import numpy as np
import matplotlib.pyplot as plt
import math

path = "eigenfunction.out"
d = np.loadtxt(path, delimiter=";")

# initialise plot
fig, ax = plt.subplots(2,1)
fig.tight_layout()

SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20

# potential 
ax[0].plot(d[:,0],d[:,1], "b", linewidth=3)
ax[0].plot(d[:, 0], d[:, 2], "r--", linewidth=3)
ax[0].plot(d[:, 0], d[:, 3]* d[-1,1]/d[-1,3], "g--", linewidth=3)
# ax[0].plot(d[:, 0], d[:, 4], "k", linewidth=3)
ax[0].plot(d[:, 0], d[:, 6] * d[-1,1]/d[-1,6], "m-.", linewidth=3)
ax[0].plot(d[:, 0], d[:, 7] * d[-1,1]/d[-1,7], "r-.", linewidth=3)
ax[0].legend([ 'Exact','Unperturbed','Basis sum','Galerkin','Augment Galerkin'], fontsize = MEDIUM_SIZE)
ax[0].tick_params(axis='both', which='major', labelsize=SMALL_SIZE)
t = ax[0].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)

ax[1].plot(d[:,0],abs(d[:,3]* d[-1,1]/d[-1,3] - d[:,1])/abs(max(d[:,1])) * 100.0, "b", linewidth=3)
ax[1].plot(d[:,0],abs(d[:,4] - d[:,1])/abs(max(d[:,1])) * 100.0, "r", linewidth=3)
ax[1].plot(d[:,0],abs(d[:,6]* d[-1,1]/d[-1,6] - d[:,1])/abs(max(d[:,1])) * 100.0, "g--", linewidth=3)
ax[1].plot(d[:,0],abs(d[:,7]* d[-1,1]/d[-1,7] - d[:,1])/abs(max(d[:,1])) * 100.0, "k", linewidth=3)
ax[1].tick_params(axis='both', which='major', labelsize=MEDIUM_SIZE)
ax[1].legend(['Standard','Augment','Galerkin','Augmented Galerkin'], fontsize=MEDIUM_SIZE)
t = ax[1].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)

plt.show()