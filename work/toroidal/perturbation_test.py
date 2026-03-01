#%%
import numpy as np
import matplotlib.pyplot as plt
import math


path = "traction_pert.out"
d = np.loadtxt(path, delimiter=";")
path2 = "eigenvectors_pert.out"
d2 = np.loadtxt(path2, delimiter=";")

# path2 = "augment_f2.out"
# d2 = np.loadtxt(path2, delimiter=";")
# d1 = np.loadtxt(path1, delimiter=";")

# initialise plot
# maxval = 4
nx = d.shape[0]
numc = nx - 5
print(nx)

# len_d = d.shape[0]

# fig, ax = plt.subplots(2,2)

# for i in range(0,2):
#     ax[0,0].plot(d[:,0],d[:,5 * i + 1],'-',linewidth=3)
#     ax[0,0].plot(d[:,0],d[:,5* i + 2],'--',linewidth=3)
#     # ax[0].plot(d[:,0],d[:,4 * i + 3],'-.',linewidth=3)
#     # ax[1].plot(d[nx-20:nx-1,0],d[nx-20:nx-1,2*i+2] - d[nx-20:nx-1,2*i+1],linewidth=3)
#     ax[1,0].plot(d[:,0],d[:,5*i+2] - d[:,5*i+1],'-',linewidth=3)
#     # ax[1].plot(d[:,0],d[:,4*i+3] - d[:,4*i+1],'--',linewidth=3)

# for i in range(0,2):
#     ax[0,1].plot(d[:,0],d[:,5 * i + 3],'-',linewidth=3)
#     ax[0,1].plot(d[:,0],d[:,5* i + 4]* np.sign(d[numc,5*i+4]/d[numc,5*i+3]),'--',linewidth=3)
#     ax[0,1].plot(d[:,0],d[:,5* i + 5] * np.sign(d[numc,5*i+5]/d[numc,5*i+3]),'-.',linewidth=3)
#     # ax[0].plot(d[:,0],d[:,4 * i + 3],'-.',linewidth=3)
#     # ax[1].plot(d[nx-20:nx-1,0],d[nx-20:nx-1,2*i+2] - d[nx-20:nx-1,2*i+1],linewidth=3)
#     ax[1,1].plot(d[:,0],d[:,5*i+3] - d[:,5*i+4]* np.sign(d[numc,5*i+4]/d[numc,5*i+3]),'-',linewidth=3)
#     ax[1,1].plot(d[:,0],d[:,5*i+3] - d[:,5*i+5]* np.sign(d[numc,5*i+5]/d[numc,5*i+3]),'--',linewidth=3)

fig, ax = plt.subplots(1,2)


for i in range(0,2):
    ax[0].plot(d[:,0],d[:,5 * i + 3],'-',linewidth=3)
    ax[0].plot(d[:,0],d[:,5* i + 4]* np.sign(d[numc,5*i+4]/d[numc,5*i+3]),'--',linewidth=3)
    ax[0].plot(d[:,0],d[:,5* i + 5] * np.sign(d[numc,5*i+5]/d[numc,5*i+3]),'-.',linewidth=3)
    # ax[0].plot(d[:,0],d[:,4 * i + 3],'-.',linewidth=3)
    # ax[1].plot(d[nx-20:nx-1,0],d[nx-20:nx-1,2*i+2] - d[nx-20:nx-1,2*i+1],linewidth=3)
    ax[1].plot(d[:,0],d[:,5*i+3] - d[:,5*i+4]* np.sign(d[numc,5*i+4]/d[numc,5*i+3]),'-',linewidth=3)
    ax[1].plot(d[:,0],d[:,5*i+3] - d[:,5*i+5]* np.sign(d[numc,5*i+5]/d[numc,5*i+3]),'--',linewidth=3)

# ax[0].set_ylabel('Traction',fontsize=20)
# ax[1].set_ylabel('Traction difference',fontsize=20)
# ax[1].set_xlabel('Radius', fontsize = 20)

# ax[0].plot(d2[:,0],d2[:,1],'b',linewidth=3)
# ax[0].plot(d2[:,0],d2[:,2],'r',linewidth=3)
# ax[0].plot(d2[:,0],d2[:,3] * d2[-1,2]/d2[-1,3],'k--',linewidth=3)
# ax.set_yscale('log')

plt.show()