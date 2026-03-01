#%%
import numpy as np
import matplotlib.pyplot as plt
import math

# d = np.loadtxt("Lagrange.out", delimiter=";")
# path_to_folder = "Bench1"
path1 = "prempert.out"
path2 = "prempert2.out"
path3 = "premaug.out"
# path2 = path_to_folder + "/IntegralSolution.out"
# path3 = path_to_folder + "/MatrixSolution.out"
# d = np.loadtxt("cleanbench1.out", delimiter=";")
d0 = np.loadtxt(path1, delimiter=";")
d1 = np.loadtxt(path2, delimiter=";")
d2 = np.loadtxt(path3, delimiter=";")
# d1 = np.loadtxt(path2, delimiter=";")
# d2 = np.loadtxt(path3, delimiter=";")
# d2 = np.loadtxt('time_test.STAT1.Z')

# initialise plot
fig, ax = plt.subplots(2,1)
# plt.rcParams.update({'font.size': 18})
fig.tight_layout()

SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20

# row, col = d1.shape
# l = int(math.sqrt((col - 1)/2)-1)
# print("Number of cols: ", col)
# print("Check: ", math.sqrt((col - 1)/2))
# print("l: ", l)

# potential 
ax[0].plot(d0[:,0],d0[:,1], "b", linewidth=3)
ax[0].plot(d0[:, 0], d0[:, 2], "r", linewidth=3)
ax[0].plot(d0[:, 0], d0[:, 3], "g-.", linewidth=3)
# ax[0].plot(d0[:, 0], d0[:, 4], "m", linewidth=3)

# ax[0].plot(d1[:, 0], d1[:, 2], "c-.", linewidth=3)
# ax[0].plot(d1[:, 0], d1[:, 3], "k-.", linewidth=3)
# ax[0].plot(d0[:,0], d0[:,2] - d0[:,1],'k',linewidth=3)
# ax[0].plot(d0[:, 0], d0[:, 5], "k", linewidth=3)

ax[0].plot(d2[:,0],d2[:,3],'r--',linewidth=3)
# ax[0].plot(d2[:, 0], d2[:, 1], "g-.", linewidth=3)
# ax[0].plot(d0[:, 0], d0[:, 1], "ko")
# ax[0].legend([ 'Exact perturbed','Basis sum','PREM','Lagrange','Rad mesh pert', 'Rad mesh sum','199th mode'], fontsize = MEDIUM_SIZE)
ax[0].legend([ 'Exact perturbed','Basis sum','PREM','augment'], fontsize = MEDIUM_SIZE)
ax[0].tick_params(axis='both', which='major', labelsize=SMALL_SIZE)
t = ax[0].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)
ax[0].set_ylabel("U", fontsize = BIGGER_SIZE)

# error in potential
ax[1].plot(d0[:,0], abs(d0[:,2] - d0[:,1])/abs(max(d0[:,1])) * 100.0, "b", linewidth=3)
ax[1].plot(d2[:,0], abs(d2[:,2] - d2[:,3])/abs(max(d0[:,2])) * 100.0, "r", linewidth=3)
ax[1].plot(d0[:,0], abs(d0[:,3] - d0[:,1])/abs(max(d0[:,1])) * 100.0, "g", linewidth=3)
# ax[1].plot(d0[:,0], abs(d0[:,4] - d0[:,1])/abs(max(d0[:,1])) * 100.0, "r", linewidth=3)
# ax[1].plot(d1[:,0], abs(d1[:,2] - d1[:,3])/abs(max(d1[:,2])) * 100.0, "k", linewidth=3)
# ax[1].plot(d0[:,0], abs(d0[:,4] - d0[:,1])/abs(max(d0[:,1])) * 100.0, "r", linewidth=3)
# ax[1].plot(d0[:,0], abs(d0[:,4] - d0[:,1])/abs(max(d0[:,1])) * 100.0, "r", linewidth=3)
# ax[1].plot(d0[:,0], (abs(d0[:,2] - d0[:,1])/abs(max(d0[:,1])) - abs(d0[:,4] - d0[:,1])/abs(max(d0[:,1])) ) * 100.0, "g", linewidth=3)
# ax[1].plot(d0[:,0], abs(d2[:,1] - d0[:,1])/abs(max(d0[:,1])) * 100.0, "r", linewidth=3)
# ax[1].legend(['Relative error (%)','Lagrange','Rad mesh'],fontsize = MEDIUM_SIZE)
ax[1].tick_params(axis='both', which='major', labelsize=MEDIUM_SIZE)
t = ax[1].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)
ax[1].set_xlabel("Radius (scaled to Earth's)", fontsize = BIGGER_SIZE)


# plt.ylim(-0.00001, 0.00001)
# manager = plt.get_current_fig_manager()
# manager.full_screen_toggle()
#figure = plt.gcf()  
#figure.set_size_inches(9, 12) 
#plt.savefig("../../../Notes/New_coupling/figures/Benchmark_1.pdf",format = "pdf", bbox_inches='tight', pad_inches = 0)
plt.show()

