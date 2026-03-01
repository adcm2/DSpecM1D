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

# path1 = "../../../YSpec/build/examples/yspec.out_spectra.1"
# d1 = np.loadtxt(path1, delimiter=" ")
path1 = "groundresponse_l_transverse.out"
d1 = np.loadtxt(path1, delimiter=";")

# path = "yspec.out_spectra.1"   # adjust as needed
data1 = read_spectra_to_array(path1)

# initialise plot
fig, ax = plt.subplots(3,3)
plt.subplots_adjust(wspace=0, hspace=0)
SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20    


ax[0,0].plot(d1[:,0],d1[:,12], "b", linewidth=3)
ax[0,0].plot(d1[:,0],d1[:,3], "g--", linewidth=3)
ax[0,0].plot(d1[:,0],d1[:,21], "k--", linewidth=3)
ax[0,0].legend(['YSpec','Transverse','MINEOS'], loc='upper right', fontsize=BIGGER_SIZE)

ax[1,0].plot(d1[:,0],np.abs(d1[:,3]- d1[:,12])/max(np.abs(d1[:,12]))*100, "k", linewidth=3)
ax[2,0].plot(d1[:,0],np.abs(d1[:,3] - d1[:,21])/max(np.abs(d1[:,3]))*100, "g--", linewidth=3)

ax[1,0].set_ylabel("Vs YSpec (%)", fontsize=MEDIUM_SIZE)
ax[2,0].set_ylabel("Vs MINEOS (%)", fontsize=MEDIUM_SIZE)



ax[0,1].plot(d1[:,0],d1[:,15], "b", linewidth=3)
ax[0,1].plot(d1[:,0],d1[:,6], "g--", linewidth=3)
ax[0,1].plot(d1[:,0],d1[:,24], "k--", linewidth=3)

# ax[1,1].plot(d1[:,0],np.abs(d1[:,15] - d1[:,6])/max(d1[:,4])*100, "b", linewidth=3)
ax[1,1].plot(d1[:,0],np.abs(d1[:,6] - d1[:,15])/max(d1[:,15])*100, "k", linewidth=3)
ax[2,1].plot(d1[:,0],np.abs(d1[:,6] - d1[:,24])/max(d1[:,24])*100, "g--", linewidth=3)
# ax[2,1].plot(d1[:,0],np.abs(d1[:,15] - d1[:,24])/max(d1[:,24])*100, "g--", linewidth=3)
# ax[1,1].legend(['Mine vs MINEOS','YSpec vs MINEOS'], loc='upper right', fontsize=BIGGER_SIZE)
# ax[1,1].set_ylabel("Percent  (%)", fontsize=MEDIUM_SIZE)

ax[0,2].plot(d1[:,0],d1[:,18], "b", linewidth=3)
ax[0,2].plot(d1[:,0],d1[:,9], "g--", linewidth=3)
ax[0,2].plot(d1[:,0],d1[:,27], "k--", linewidth=3)
# ax[1,2].plot(d1[:,0],np.abs(d1[:,18] - d1[:,9])/max(d1[:,9])*100, "b", linewidth=3)
ax[1,2].plot(d1[:,0],np.abs(d1[:,9] - d1[:,18])/max(d1[:,18])*100, "k", linewidth=3)
ax[2,2].plot(d1[:,0],np.abs(d1[:,9] - d1[:,27])/max(d1[:,27])*100, "g--", linewidth=3)
# ax[2,2].plot(d1[:,0],np.abs(d1[:,18] - d1[:,27])/max(d1[:,27])*100, "g--", linewidth=3)
# ax[1,2].legend(['Mine vs MINEOS','YSpec vs MINEOS'], loc='upper right', fontsize=BIGGER_SIZE)
# ax[1,2].set_ylabel("Percent Error (%)", fontsize=MEDIUM_SIZE)

ax[0,0].set_title("Vertical Ground Motion")
ax[0,1].set_title("North-South Ground Motion")
ax[0,2].set_title("East-West Ground Motion")
plt.show()