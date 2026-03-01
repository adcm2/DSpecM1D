#%%
import numpy as np
import matplotlib.pyplot as plt
import math

# path = "benchmark.out"
# d = np.loadtxt(path, delimiter=";")

fig, ax = plt.subplots(1,1)

path0 = "mineos_bench.out"
for i in range(1,11):
    path2 = path0 + "_" + str(i)
    d1 = np.loadtxt(path2,delimiter=";")
    ax.plot(d1[:,0],d1[:,1],marker="o",markersize=8,linestyle='-',linewidth=3)

leginfo = list()
for i in range(1,11):
    strl = 'l = ' + str(i)
    leginfo.append(strl)

ax.set_yscale('log')
ax.set_xlabel('n', fontsize = 20)
ax.set_ylabel('Relative difference', fontsize = 20)
ax.legend(leginfo, fontsize = 15, loc ='upper right')
# ax[0].legend([ 'Exact','Matrix-free pseudospectral','Nodes'], fontsize = MEDIUM_SIZE)
ax.tick_params(axis='both', which='major', labelsize=20)
# t = ax[0].yaxis.get_offset_text()
# t.set_size(MEDIUM_SIZE)
# ax[0].set_ylabel("Potential", fontsize = BIGGER_SIZE)
plt.show()