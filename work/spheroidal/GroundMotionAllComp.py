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

path = "groundresponse_l.out"
d = np.loadtxt(path, delimiter=";")
path1 = "groundresponse_l_transverse.out"
d1 = np.loadtxt(path1, delimiter=";")
path2 = "groundresponse_raw.out"
d2 = np.loadtxt(path2, delimiter=";")


# initialise plot
maxval = 4
fig, ax = plt.subplots(1,3)
plt.subplots_adjust(wspace=0, hspace=0)
SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20    

ax[0].plot(d[:,0],d[:,12], "b", linewidth=3)
ax[0].plot(d[:,0],d[:,3], "g--", linewidth=3)
ax[0].plot(d1[:,0],d1[:,3], "r--", linewidth=3)
ax[0].legend(['YSpec','Isotropic','Transverse'], loc='upper right', fontsize=BIGGER_SIZE)

ax[1].plot(d[:,0],d[:,15], "b", linewidth=3)
ax[1].plot(d[:,0],d[:,6], "g--", linewidth=3)
ax[1].plot(d1[:,0],d1[:,6], "r--", linewidth=3)
ax[1].legend(['YSpec','Isotropic','Transverse'], loc='upper right', fontsize=BIGGER_SIZE)

ax[2].plot(d[:,0],d[:,18], "b", linewidth=3)
ax[2].plot(d[:,0],d[:,9], "g--", linewidth=3)
ax[2].plot(d1[:,0],d1[:,9], "r--", linewidth=3)
ax[2].legend([ 'YSpec','Isotropic','Transverse'], loc='upper right', fontsize=BIGGER_SIZE)



ax[0].set_title('Z amplitude', fontsize=BIGGER_SIZE)
ax[1].set_title('Theta amplitude', fontsize=BIGGER_SIZE)
ax[2].set_title('Phi amplitude', fontsize=BIGGER_SIZE)

ax[0].set_xlabel('Frequency (mHz)', fontsize=BIGGER_SIZE)

plt.show()