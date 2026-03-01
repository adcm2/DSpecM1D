#%%
import numpy as np
import matplotlib.pyplot as plt
import math

path = "functiontozero.out"
d = np.loadtxt(path, delimiter=";")

fig, ax = plt.subplots(1,1)

ax.plot(d[:,0],d[:,1], "b", linewidth=3)
plt.show()