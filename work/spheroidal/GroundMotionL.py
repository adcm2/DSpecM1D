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

path_mine = "groundresponse_raw_transverse.out"
d_mine = np.loadtxt(path_mine, delimiter=";")

path1 = "../../../YSpec/output/yspec.out_spectra.1"
data1 = read_spectra_to_array(path1)

path_filt = "../../../YSpec/output/yspec.out_spectra_filt.1"
data_filt = read_spectra_to_array(path_filt)

path_mfilt = "check_filt.out"
d_mfilt = np.loadtxt(path_mfilt, delimiter=";")

path_yspec_t = "../../../YSpec/output/yspec.out.1"
data_yspec_t = read_spectra_to_array(path_yspec_t)

path_t = "check_time.out"
data_t = np.loadtxt(path_t, delimiter=";")

# path = "yspec.out_spectra.1"   # adjust as needed

# print("shape:", data1.shape)
# print("first row:", data1[0])
# print("first row:", data1[1])

# initialise plot
maxval = 4
fig, ax = plt.subplots(2,2)
plt.subplots_adjust(wspace=0, hspace=0)
SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20    

print(data1.shape)
print(d_mine.shape)
print(data1[0:5,0])
print(d_mine[0:5,0])
# ax[0,0].plot(data1[:,0],data1[:,1], "b", linewidth=3)
# ax[0,0].plot(d_mfilt[:,0],d_mfilt[:,1], "r--", linewidth=3)
# ax[0,0].legend(["YSpec unfiltered","My code unfiltered"], fontsize=SMALL_SIZE)

# ax[0,0].plot(data1[:,0],data1[:,1], "b", linewidth=3)
# ax[0,0].plot(d_mine[:,0],d_mine[:,1], "r--", linewidth=3)
# ax[0,0].plot()
# ax[1,0].plot(data1[:,0],data1[:,2], "b", linewidth=3)
# ax[1,0].plot(d_mine[:,0],d_mine[:,2], "r--", linewidth=3)

ax[0,0].plot(data_filt[:,0],data_filt[:,1], "b", linewidth=3)
ax[0,0].plot(d_mfilt[:,0],d_mfilt[:,3], "r--", linewidth=3)

ax[1,0].plot(data_filt[:,0],100 * abs(data_filt[:,1] - d_mfilt[:,1])/max(abs(d_mfilt[:,1])), "b", linewidth=3)

ax[1,0].legend(["YSpec filtered","My code filtered"], fontsize=SMALL_SIZE)
ax[1,0].set_xlabel("Frequency (mHz)", fontsize=MEDIUM_SIZE)

ax[0,1].plot(data_yspec_t[:,0],data_yspec_t[:,1], "b", linewidth=3)#
ax[0,1].plot(data_t[:,0],data_t[:,1], "r--", linewidth=3)
ax[0,1].legend(["YSpec time series","My code time series"], fontsize=SMALL_SIZE)

# ax[1,1].plot(data_yspec_t[:,0],100 * (data_yspec_t[:,1] - data_t[:,1])/max(data_t[:,1]), "b", linewidth=3)

ax[1,1].set_xlabel("Time (s)", fontsize=MEDIUM_SIZE)
   
plt.show()