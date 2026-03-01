#%%
import numpy as np
import matplotlib.pyplot as plt
import math

path = "benchmark.out"
d = np.loadtxt(path, delimiter=";")

fig, ax = plt.subplots(1,1)

ax.plot(d[:,1],d[:,2], "bo", linewidth=3)
ax.plot(d[:,3],d[:,4], "co", linewidth=3)
ax.plot(d[:,5],d[:,6], "go", linewidth=3)
ax.plot(d[:,7],d[:,8], "mo", linewidth=3)
ax.plot(d[:,9],d[:,10], "ro", linewidth=3)
ax.plot(d[:,11],d[:,12], "ko", linewidth=3)

z1 = np.polyfit(np.log(d[:,1]),np.log(d[:,2]),1)
z2 = np.polyfit(np.log(d[:,3]),np.log(d[:,4]),1)
z3 = np.polyfit(np.log(d[:,5]),np.log(d[:,6]),1)
z4 = np.polyfit(np.log(d[:,7]),np.log(d[:,8]),1)
z5 = np.polyfit(np.log(d[:,9]),np.log(d[:,10]),1)
z6 = np.polyfit(np.log(d[:,11]),np.log(d[:,12]),1)

# z3 = np.polyfit(np.log(d[0:4,0]),np.log(d[0:4,3]),1)
# z4 = np.polyfit(np.log(d[0:3,0]),np.log(d[0:3,4]),1)
# z5 = np.polyfit(np.log(d[0:2,0]),np.log(d[0:2,5]),1)
# z6 = np.polyfit(np.log(d[0:2,0]),np.log(d[0:2,6]),1)

print(z1)
print(z2)
print(z3)
print(z4)
print(z5)
print(z6)

ax.plot(d[:,1], np.exp(z1[1]) * np.power(d[:,1],z1[0]),"b",linewidth=2)
ax.plot(d[:,3], np.exp(z2[1]) * np.power(d[:,3],z2[0]),"c",linewidth=2)
ax.plot(d[:,5], np.exp(z3[1]) * np.power(d[:,5],z3[0]),"g",linewidth=2)
ax.plot(d[:,7], np.exp(z4[1]) * np.power(d[:,7],z4[0]),"m",linewidth=2)
ax.plot(d[:,9], np.exp(z5[1]) * np.power(d[:,9],z5[0]),"r",linewidth=2)
ax.plot(d[:,11], np.exp(z6[1]) * np.power(d[:,11],z6[0]),"k",linewidth=2)
# ax.plot(d[0:6,0], np.exp(z2[1]) * np.power(d[0:6,0],z2[0]),"c",linewidth=2)
# ax.plot(d[0:4,0], np.exp(z3[1]) * np.power(d[0:4,0],z3[0]),"g",linewidth=2)
# ax.plot(d[0:3,0], np.exp(z4[1]) * np.power(d[0:3,0],z4[0]),"m",linewidth=2)
# ax.plot(d[0:2,0], np.exp(z5[1]) * np.power(d[0:2,0],z5[0]),"r",linewidth=2)
# ax.plot(d[0:2,0], np.exp(z6[1]) * np.power(d[0:2,0],z6[0]),"k",linewidth=2)
# ax.plot(d[:,0],np.power(d[:,1],6) * d[0,1]/np.power(d[0,0],6),"k--", linewidth = 2)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('h (arb. units)', fontsize = 20)
ax.set_ylabel('Relative error', fontsize = 20)
ax.legend(['p = 3', 'p = 4','p = 5', 'p = 6','p = 7', 'p = 8'], fontsize = 20)
# ax[0].legend([ 'Exact','Matrix-free pseudospectral','Nodes'], fontsize = MEDIUM_SIZE)
ax.tick_params(axis='both', which='major', labelsize=20)
# t = ax[0].yaxis.get_offset_text()
# t.set_size(MEDIUM_SIZE)
# ax[0].set_ylabel("Potential", fontsize = BIGGER_SIZE)
plt.show()