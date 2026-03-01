#%%
import numpy as np
import matplotlib.pyplot as plt
import math

path = "groundresponse.out"
d = np.loadtxt(path, delimiter=";")

# initialise plot
maxval = 4
fig, ax = plt.subplots(2,2)
plt.subplots_adjust(wspace=0, hspace=0)
SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20    
# ax[0].plot(d[:,0],d[:,3], "b", linewidth=3)
# ax[0].plot(d[:,0],d[:,6], "r", linewidth=3)
# ax[0].plot(d[:,0],d[:,9], "m--", linewidth=3)
# ax[0].plot(d[:,0],d[:,10], "y--", linewidth=3)
ax[0,0].plot(d[:,0],d[:,1], "g", linewidth=3)
ax[0,0].plot(d[:,0],d[:,11], "c--", linewidth=3)
ax[0,0].plot(d[:,0],d[:,2], linewidth=3)
ax[0,0].plot(d[:,0],d[:,13], "--", linewidth=3)
ax[0,0].plot(d[:,0],d[:,15], "k--", linewidth=3)

ax[0,1].plot(d[:,0],d[:,4], "g", linewidth=3)
ax[0,1].plot(d[:,0],d[:,12], "c--", linewidth=3)
ax[0,1].plot(d[:,0],d[:,5], linewidth=3)
ax[0,1].plot(d[:,0],d[:,14], "--", linewidth=3)
# ax[0].plot(d[:,0],d[:,13], "k-.", linewidth=3)
# ax[0].plot(d[:,0],d[:,14], "k", linewidth=3)
# ax[1].plot(d[:,0],abs(d[:,3] - d[:,7])/max(d[:,3])*100, "b", linewidth=3)
# ax[1].plot(d[:,0],abs(d[:,6] - d[:,8])/max(d[:,6])*100, "r", linewidth=3)
# ax[1].plot(d[:,0],abs(d[:,3] - d[:,9])/max(d[:,3])*100, "g", linewidth=3)
# ax[1].plot(d[:,0],abs(d[:,6] - d[:,10])/max(d[:,6])*100, "m", linewidth=3)

ax[1,0].plot(d[:,0],abs(d[:,1] - d[:,11])/max(d[:,1])*100, "c", linewidth=3)
ax[1,0].plot(d[:,0],abs(d[:,2] - d[:,13])/max(d[:,2])*100, "k", linewidth=3)

ax[1,1].plot(d[:,0],abs(d[:,4] - d[:,12])/max(d[:,4])*100, "c", linewidth=3)
ax[1,1].plot(d[:,0],abs(d[:,5] - d[:,14])/max(d[:,5])*100, "k", linewidth=3)
# ax[1].plot

ax[1,0].set_xlabel('Frequency (mHz)', fontsize=BIGGER_SIZE)
ax[1,0].set_ylabel('Relative error (%)', fontsize=BIGGER_SIZE)
ax[0,0].set_ylabel('Amplitude', fontsize=BIGGER_SIZE)
ax[0,0].set_xlabel('Frequency (mHz)', fontsize=BIGGER_SIZE)
# ax[0].plot( d2[:,0],d2[:,i] * np.sign(d[-1,i]/d2[-1,i]),"r--", linewidth=3)
    # ax[0][0,i-1].plot(d3[:,0],d3[:,i], "k", linewidth=3)
    # ax[0][0,i-1].plot(d4[:,0],d4[:,i]* np.sign(d4[-1,i]/d3[-1,i]), "g", linewidth=3)
    # ax[0][0,i-1].plot(d5[:,0],d5[:,i], "m", linewidth=3)
    # ax[0][0,i-1].plot(d6[:,0],d6[:,i], "y", linewidth=3)
    
    # ax[0][1,i-1].plot(d2[:,0],abs(d2[:,i] * np.sign(d[-1,i]/d2[-1,i]) - d[:,i])/max[0](d[:,i]) * 100,'b')
    # ax[0][1,i-1].plot(d4[:,0],abs(d4[:,i] * np.sign(d4[-1,i]/d3[-1,i]) - d3[:,i])/max[0](d4[:,i]) * 100,'k')
# ax[0].legend(['Gen eig sol', 'Mineos'])


# ax[0].legend(['Theta motion abs', 'Phi motion abs', 'T2', 'P2', 'NMCT', 'NMCP'], loc='upper right', fontsize=BIGGER_SIZE)      
# ax[1].legend(['Theta error', 'Phi error', 'Theta error NMCT', 'Phi error NMCP'], loc='upper right', fontsize=BIGGER_SIZE)

ax[1,0].legend(['Real error', 'Imag error'], loc='upper right', fontsize=BIGGER_SIZE)
ax[1,1].legend(['Real error', 'Imag error'], loc='upper right', fontsize=BIGGER_SIZE)

# ax[0].set_xlabel('Frequency (mHz)', fontsize=BIGGER_SIZE)
# ax[0].set_ylabel('Amplitude', fontsize=BIGGER_SIZE)      
plt.show()