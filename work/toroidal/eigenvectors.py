#%%
import numpy as np
import matplotlib.pyplot as plt
import math


path = "eigenvectors_pert.out"
d = np.loadtxt(path, delimiter=";")

fig, ax = plt.subplots(2,2)

nc = d.shape[0]-1

for i in range(0,2):
    ax[0,i].plot(d[:,0],d[:,4* i + 1],'b',linewidth=3)
    ax[0,i].plot(d[:,0],d[:,4* i + 1],'bo',linewidth=3)
    ax[0,i].plot(d[:,0],d[:,4* i + 2],'r',linewidth=3)
    ax[0,i].plot(d[:,0],d[:,4* i + 3]*np.sign(d[nc,4*i+3]/d[nc,4*i+1]),'g--',linewidth=3)
    ax[0,i].plot(d[:,0],d[:,4* i + 4]*np.sign(d[nc,4*i+4]/d[nc,4*i+1]),'k-.',linewidth=3)
    ax[0,i].legend(['Original','Perturbed','Unaugmented','Augmented'])

    # ax[1,i].plot(d[:,0],abs(d[:,4* i + 1]-d[:,4* i + 2]),'r',linewidth=3)
    ax[1,i].plot(d[:,0],abs(d[:,4* i + 3]*np.sign(d[nc,4*i+3]/d[nc,4*i+1])-d[:,4*i+2]),'g-',linewidth=3)
    ax[1,i].plot(d[:,0],abs(d[:,4* i + 4]*np.sign(d[nc,4*i+4]/d[nc,4*i+1])-d[:,4*i+2]),'k-',linewidth=3)
    ax[1,i].legend(['Unaugmented','Augmented'])


plt.show()