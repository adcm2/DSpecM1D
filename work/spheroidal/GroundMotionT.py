#%%
import numpy as np
import matplotlib.pyplot as plt
import math
import re

_float_re = re.compile(r'[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][+-]?\d+)?')

def read_spectra_to_array(path, pad=np.nan):
    """
    Read numeric columns from a text spectra file into a 2D numpy array.
    Lines with different column counts are padded with `pad` (default NaN).
    """
    rows = []
    with open(path, 'r') as fh:
        for line in fh:
            nums = _float_re.findall(line)
            if not nums:
                continue
            rows.append([float(x) for x in nums])
    if not rows:
        return np.empty((0, 0), dtype=float)
    maxcols = max(len(r) for r in rows)
    arr = np.full((len(rows), maxcols), pad, dtype=float)
    for i, r in enumerate(rows):
        arr[i, :len(r)] = r
    return arr

# path = "groundresponse_l.out"
# d = np.loadtxt(path, delimiter=";")
# path2 = "groundresponse_raw.out"
# d2 = np.loadtxt(path2, delimiter=";")

path = "groundresponse_time.out"
d = np.loadtxt(path, delimiter=";")
path1 = "groundresponse_t_yspec.out"
d2 = np.loadtxt(path1, delimiter=";")

# initialise plot
maxval = 4
fig, ax = plt.subplots(1)
plt.subplots_adjust(wspace=0, hspace=0)
SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20    

ax.plot(d[:,0]/3600,d[:,1], "b", linewidth=3)
# ax.plot(d2[:,0],d2[:,3]/max(d[:,3]), "k", linewidth=3)
# ax[0].plot(d[:,0],d[:,2], "r", linewidth=3)
# ax.plot(d[:,0],d[:,3],"k", linewidth=3)
# mag_ur = np.hypot(d2[:,1], d2[:,2])
# ax.plot(data1[:,0], data1[:,1], "g--", linewidth=2)
# ax.plot(data1[:,0], data1[:,2], "m--", linewidth=2)
ax.plot(d2[:,0]/3600, d2[:,1], "c--", linewidth=2)

# ax[1].plot(d[:,0],d[:,4], "b", linewidth=3)
# ax[1].plot(d[:,0],d[:,5], "r", linewidth=3)
# ax[1].plot(d[:,0],d[:,6],"k", linewidth=3)

# ax[2].plot(d[:,0],d[:,7], "b", linewidth=3)
# ax[2].plot(d[:,0],d[:,8], "r", linewidth=3)
# ax[2].plot(d[:,0],d[:,9],"k", linewidth=3)

ax.set_title('Z amplitude', fontsize=BIGGER_SIZE)
# ax[1].set_title('Theta amplitude', fontsize=BIGGER_SIZE)
# ax[2].set_title('Phi amplitude', fontsize=BIGGER_SIZE)

ax.set_xlabel('Time (hr)', fontsize=BIGGER_SIZE)
# ax.set_yscale('log')
# ax[1].set_xlabel('Frequency (mHz)', fontsize=BIGGER_SIZE)
# ax[2].set_xlabel('Frequency (mHz)', fontsize=BIGGER_SIZE)

# ax[0].legend(['Real Z', 'Imag Z', 'Abs Z'], loc='upper right', fontsize=BIGGER_SIZE)
# ax[1].legend(['Real Theta', 'Imag Theta', 'Abs Theta'], loc='upper right', fontsize=BIGGER_SIZE)
# ax[2].legend(['Real Phi', 'Imag Phi', 'Abs Phi'], loc='upper right', fontsize=BIGGER_SIZE)

# ax[0].set_xlabel('Frequency (mHz)', fontsize=BIGGER_SIZE)
# ax[0].set_ylabel('Amplitude', fontsize=BIGGER_SIZE)      
plt.show()