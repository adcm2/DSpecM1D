#%%
import numpy as np
import matplotlib.pyplot as plt
import math

path = "eigenfunction.out"
path1 = "eigenfunction_spec.out"
d = np.loadtxt(path, delimiter=";")
d1 = np.loadtxt(path1, delimiter=";")

# initialise plot
fig, ax = plt.subplots(2,1)
fig.tight_layout()

SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20

# potential 
ax[0].plot(d[:,0],d[:,2], "b", linewidth=3)
ax[0].plot(d1[:,0],d1[:,1] * d[-1,2]/d1[-1,1], "r--", linewidth=3)
ax[0].plot(d1[:,0],d1[:,3] * d[-1,2]/d1[-1,3], "g--", linewidth=3)
# ax[0].plot(d[:, 0], d[:, 2], "r--", linewidth=3)
# ax[0].plot(d[:, 0], d[:, 3], "g-.", linewidth=3)
# ax[0].plot(d[:, 0], d[:, 4], "k", linewidth=3)
# ax[0].plot(d[:, 0], d[:, 5], "m", linewidth=3)
ax[0].legend([ 'Exact','Spectral element','Galerkin'], fontsize = MEDIUM_SIZE)
ax[0].set_title('Comparison of unperturbed eigenfunction',fontsize=MEDIUM_SIZE)
ax[0].tick_params(axis='both', which='major', labelsize=SMALL_SIZE)
t = ax[0].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)

ax[1].plot(d[:,0],abs(d[:,2] - d1[:,1]* d[-1,2]/d1[-1,1])/abs(max(d[:,2])) * 100.0, "b", linewidth=3)
ax[1].plot(d[:,0],abs(d[:,2] - d1[:,3] * d[-1,2]/d1[-1,3])/abs(max(d[:,2])) * 100.0, "r--", linewidth=3)
# ax[1].plot(d[:,0],abs(d[:,4] - d[:,1])/abs(max(d[:,1])) * 100.0, "r", linewidth=3)
ax[1].tick_params(axis='both', which='major', labelsize=MEDIUM_SIZE)
ax[1].legend(['Relative error (%)', 'Galerkin'], fontsize=MEDIUM_SIZE)
t = ax[1].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)

plt.show()