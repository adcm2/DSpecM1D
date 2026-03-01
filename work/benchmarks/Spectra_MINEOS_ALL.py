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
# path1 = "f_domain_maxstep_100_f11_0.1_f12_0.1_f21_0.1_f22_0.1.out"
# path1 = "f_domain_maxstep_0.005000_f11_0.100000_f12_0.150000_f21_4.900000_f22_5.000000.out"
path1 = "f_domain_maxstep_0.005000_f11_0.100000_f12_0.150000_f21_48.900000_f22_49.000000.out"
d1 = np.loadtxt(path1, delimiter=";")

# path = "yspec.out_spectra.1"   # adjust as needed
data1 = read_spectra_to_array(path1)

# initialise plot
fig, ax = plt.subplots(2,3)
plt.subplots_adjust(wspace=0, hspace=0)
SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20    


ax[0,0].plot(d1[:,0],d1[:,12], "b", linewidth=3)
ax[0,0].plot(d1[:,0],d1[:,3], "g-", linewidth=3)
ax[0,0].plot(d1[:,0],d1[:,21], "k--", linewidth=3)
ax[0,0].legend(['SPEC','MINEOS'], loc='upper right', fontsize=BIGGER_SIZE)
ax[1,0].plot(d1[:,0],np.abs(d1[:,3] - d1[:,21])/max(np.abs(d1[:,3]))*100, "g--", linewidth=3)
ax[1,0].set_ylabel("Vs YSpec (%)", fontsize=MEDIUM_SIZE)



ax[0,1].plot(d1[:,0],d1[:,6], "g-", linewidth=3)
ax[0,1].plot(d1[:,0],d1[:,24], "k--", linewidth=3)
ax[1,1].plot(d1[:,0],np.abs(d1[:,4] - d1[:,22])/max(d1[:,22])*100, "g--", linewidth=3)

ax[0,2].plot(d1[:,0],d1[:,9], "g-", linewidth=3)
ax[0,2].plot(d1[:,0],d1[:,27], "k--", linewidth=3)
ax[1,2].plot(d1[:,0],np.abs(d1[:,9] - d1[:,27])/max(d1[:,27])*100, "g--", linewidth=3)

ax[0,0].set_title("Vertical Ground Motion")
ax[0,1].set_title("North-South Ground Motion")
ax[0,2].set_title("East-West Ground Motion")
plt.show()