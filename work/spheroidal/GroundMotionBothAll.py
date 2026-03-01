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
path2 = "groundresponse_raw.out"
d2 = np.loadtxt(path2, delimiter=";")

# path1 = "../../../YSpec/build/examples/yspec.out_spectra.1"
# d1 = np.loadtxt(path1, delimiter=" ")
path1 = "groundresponse_l_transverse.out"
d1 = np.loadtxt(path1, delimiter=";")

# path = "yspec.out_spectra.1"   # adjust as needed
data1 = read_spectra_to_array(path1)
print("shape:", data1.shape)
print("first row:", data1[0])
print("first row:", data1[1])

# initialise plot
maxval = 4
fig, ax = plt.subplots(2,3)
plt.subplots_adjust(wspace=0, hspace=0)
SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20    

# ax[0,0].plot(d[:,0],d[:,2], "r--", linewidth=3)
ax[0,0].plot(d1[:,0],d1[:,11], "b", linewidth=3)
ax[0,0].plot(d1[:,0],d1[:,2], "g--", linewidth=3)
ax[0,0].plot(d1[:,0],d1[:,20], "k--", linewidth=3)
ax[0,0].legend(['Transverse','YSpec','MINEOS'], loc='upper right', fontsize=BIGGER_SIZE)

# ax[1,0].plot(d[:,0],d[:,3], "r--", linewidth=3)
ax[1,0].plot(d1[:,0],d1[:,12      ] , "b", linewidth=3)
ax[1,0].plot(d1[:,0],d1[:,3], "g--", linewidth=3)

# ax[0,1].plot(d[:,0],d[:,4], "r--", linewidth=3)
ax[0,1].plot(d1[:,0],d1[:,13], "b", linewidth=3)
ax[0,1].plot(d1[:,0],d1[:,4], "g--", linewidth=3)

# ax[1,1].plot(d[:,0],d[:,5], "r--", linewidth=3)
ax[1,1].plot(d1[:,0],d1[:,14], "b", linewidth=3)
ax[1,1].plot(d1[:,0],d1[:,5], "g--", linewidth=3)

# ax[0,2].plot(d[:,0],d[:,7], "r--", linewidth=3)   
ax[0,2].plot(d1[:,0],d1[:,16], "b", linewidth=3)
ax[0,2].plot(d1[:,0],d1[:,7], "g--", linewidth=3)

# ax[1,2].plot(d[:,0],d[:,8], "r--", linewidth=3)
ax[1,2].plot(d1[:,0],d1[:,17], "b", linewidth=3)
ax[1,2].plot(d1[:,0],d1[:,8], "g--", linewidth=3)

ax[0,0].set_title("Vertical Ground Motion")
ax[0,1].set_title("North-South Ground Motion")
ax[0,2].set_title("East-West Ground Motion")
plt.show()