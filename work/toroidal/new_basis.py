#%%
import numpy as np
import matplotlib.pyplot as plt
import math

# path = "efunction.out"
# d = np.loadtxt(path, delimiter=";")
# path2= "efunction_mineos.out"
# d2 = np.loadtxt(path2, delimiter=";")
# path3 = "efunctiond.out"
# d3 = np.loadtxt(path3, delimiter=";")
# path4 = "efunctiond2.out"
# d4 = np.loadtxt(path4, delimiter=";")
# path5 = "traction.out"
# d5 = np.loadtxt(path5, delimiter=";")
path6 = "augbasis.out"
d6 = np.loadtxt(path6, delimiter=";")

# path2 = "augment_f2.out"
# d2 = np.loadtxt(path2, delimiter=";")
# d1 = np.loadtxt(path1, delimiter=";")

# initialise plot
maxval = 4
fig, ax = plt.subplots(1,maxval)
# fig2, ax2 = plt.subplots(1,1)
# fig.tight_layout()

SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20

# ax = [plt.subplot(2,2,i+1) for i in range(4)]

# for a in ax:
#     a.set_xticklabels([])
#     a.set_yticklabels([])

plt.subplots_adjust(wspace=0, hspace=0)

# print(d.shape[0])
# print(d.shape[1])
nc = d6.shape[1]
print("nc is: ", nc)


# potential 
# ax2[0].plot(d[:,0],d[:,19], "b", linewidth=3)

# if l = 1
# for i in range (1,maxval + 1):
#     ax[0,i-1].plot(d[:,0],d[:,i+1], "b", linewidth=3)
#     ax[0,i-1].plot( d2[:,0],d2[:,i] * d[-1,i+1]/d2[-1,i],"r--", linewidth=3)
#     ax[1,i-1].plot(d2[:,0],abs(d2[:,i] * d[-1,i+1]/d2[-1,i] - d[:,i+1])/max(d[:,i+1]) * 100,'b')

for i in range (0,maxval):
    # ax[i].plot(d6[:,0],d6[:,3*i+1], "b", linewidth=3)
    # ax[i].plot(d6[:,0],d6[:,3*i+2], "r", linewidth=3)
    ax[i].plot(d6[:,0],d6[:,3*i+3], "g", linewidth=3)
    # ax[0,i-1].plot( d2[:,0],d2[:,i] * np.sign(d[-1,i]/d2[-1,i]),"r--", linewidth=3)
    # ax[0,i-1].plot(d3[:,0],d3[:,i], "k", linewidth=3)
    # ax[0,i-1].plot(d4[:,0],d4[:,i]* np.sign(d4[-1,i]/d3[-1,i]), "g", linewidth=3)
    # ax[0,i-1].plot(d5[:,0],d5[:,i], "m", linewidth=3)
    # ax[0,i-1].plot(d6[:,0],d6[:,i], "y", linewidth=3)
    
    # ax[1,i-1].plot(d2[:,0],abs(d2[:,i] * np.sign(d[-1,i]/d2[-1,i]) - d[:,i])/max(d[:,i]) * 100,'b')
    # ax[1,i-1].plot(d4[:,0],abs(d4[:,i] * np.sign(d4[-1,i]/d3[-1,i]) - d3[:,i])/max(d4[:,i]) * 100,'k')
# ax[0,0].legend(['Gen eig sol', 'Mineos'])
# ax[1,0].legend(['Rel diff (%)'])

# ax2.plot(d[:,0],d[:,1], "b", linewidth=3)
# ax2.plot( d2[:,0],d2[:,1] * np.sign(d[-1,1]/d2[-1,1]),"r--", linewidth=3)
# ax2.plot(d3[:,0],d3[:,1], "k", linewidth=3)
# ax2.plot(d4[:,0],d4[:,1], "g", linewidth=3)

# ax[0].plot(d2[:,1]*8,d2[:,0],"r",linewidth=3)
# for i in range(0,int(d2.shape[1]-1)):
#     # ncol = nc - int(2 * i) - 2
#     print(i)
#     ax[i].plot(d2[:,0], d2[:,i+1],"b", linewidth=3)

# ax[0].tick_params(axis='both', which='major', labelsize=SMALL_SIZE)
# t = ax[0].yaxis.get_offset_text()
# t.set_size(MEDIUM_SIZE)




    # ax[0].plot(d[:,2], d[:,0],"r", linewidth=3)
# ax[0].plot(d[:,0],d[:,1], "b", linewidth=3)
# ax[0].plot(d2[:,0],d2[:,1], "r", linewidth=3)
# ax[0].plot(d2[:,0],d2[:,2], "g", linewidth=3)
# ax[0].plot(d2[:,0],d2[:,3], "k", linewidth=3)
# ax[0].set_aspect('equal')
# ax[0].plot(d1[:,0],d1[:,1] * d[-1,2]/d1[-1,1], "r--", linewidth=3)
# ax[0].plot(d1[:,0],d1[:,3] * d[-1,2]/d1[-1,3], "g--", linewidth=3)
# ax[0].plot(d[:, 0], d[:, 2], "r--", linewidth=3)
# ax[0].plot(d[:, 0], d[:, 3], "g-.", linewidth=3)
# ax[0].plot(d[:, 0], d[:, 4], "k", linewidth=3)
# ax[0].plot(d[:, 0], d[:, 5], "m", linewidth=3)
# ax[0].legend([ 'Exact','Spectral element','Galerkin'], fontsize = MEDIUM_SIZE)
# ax[0].set_title('Comparison of unperturbed eigenfunction',fontsize=MEDIUM_SIZE)

# ax[1].plot(d[:,0],abs(d[:,2] - d1[:,1]* d[-1,2]/d1[-1,1])/abs(max(d[:,2])) * 100.0, "b", linewidth=3)
# ax[1].plot(d[:,0],abs(d[:,2] - d1[:,3] * d[-1,2]/d1[-1,3])/abs(max(d[:,2])) * 100.0, "r--", linewidth=3)
# # ax[1].plot(d[:,0],abs(d[:,4] - d[:,1])/abs(max(d[:,1])) * 100.0, "r", linewidth=3)
# ax[1].tick_params(axis='both', which='major', labelsize=MEDIUM_SIZE)
# ax[1].legend(['Relative error (%)', 'Galerkin'], fontsize=MEDIUM_SIZE)
# t = ax[1].yaxis.get_offset_text()
# t.set_size(MEDIUM_SIZE)

plt.show()