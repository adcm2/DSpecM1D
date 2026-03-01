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
# d1 = np.loadtxt(path2, delimiter=";")

# path = "groundresponse_time.out"
# d = np.loadtxt(path, delimiter=";")
path1 = "groundresponse_t_yspec_transverse.out"
d1 = np.loadtxt(path1, delimiter=";")

path2 = "groundresponse_time_transverse.out"
d2 = np.loadtxt(path2, delimiter=";")

# path_mineos_E = "/home/adcm2/space/mineos-1.0.2/DEMO/MYEX/Syndat_ASC_NOHEADER/Syndat.2000014:23:37:10.TLY.LHE.ASC"
# d_mineos_E = np.loadtxt(path_mineos_E)

# path_mineos_N = "/home/adcm2/space/mineos-1.0.2/DEMO/MYEX/Syndat_ASC_NOHEADER/Syndat.2000014:23:37:10.TLY.LHN.ASC"
# d_mineos_N = np.loadtxt(path_mineos_N)

# path_mineos_Z = "/home/adcm2/space/mineos-1.0.2/DEMO/MYEX/Syndat_ASC_NOHEADER/Syndat.2000014:23:37:10.TLY.LHZ.ASC"
# d_mineos_Z = np.loadtxt(path_mineos_Z)

path_mf = "groundresponse_t_mineos_transverse.out"
path_traw = "groundresponse_t_yspec_transverse.out"
dmf = np.loadtxt(path_mf, delimiter=";")
dtraw = np.loadtxt(path_traw, delimiter=";")

# len1 = d.shape[0]
lenl = d1.shape[0]
# lenl = min(len1, len2)

# initialise plot
maxval = 4
fig, ax = plt.subplots(3,3)
plt.subplots_adjust(wspace=0, hspace=0)
SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20    

# ax[0].plot(d[:,0]/3600,d[:,1], "b", linewidth=3)
# ax[0].plot(d1[:,0]/3600, d1[:,1], "c--", linewidth=2)
######################################################
# ax[0,0].plot(d[:,0]/3600,d[:,1], "g--", linewidth=3)
# ax[0,0].plot(d1[:,0]/3600, d1[:,1], "b", linewidth=2)
# ax[0,0].plot(d2[:,0]/3600, d2[:,1], "r--", linewidth=2)


# ax[0,0].plot(d_mineos_Z[:,0]/3600, d_mineos_Z[:,1] * 1e-9, "k", linewidth=2)
ax[0,0].plot(dmf[:,0], dmf[:,7], "b", linewidth=2)
ax[0,0].plot(dmf[:,0], dmf[:,4], "r--", linewidth=2)
ax[0,0].plot(dmf[:,0], dmf[:,1], "g--", linewidth=2)
ax[0,0].legend(['YSpec','Transverse','MINEOS'], loc='upper right', fontsize=BIGGER_SIZE)
# ax[0,0].set_xlim(400,1300)

ax[1,0].plot(dmf[:,0],(dmf[:,7]-dmf[:,4])/max(dmf[:,4]) * 100, "k", linewidth=2)
ax[2,0].plot(dmf[:,0],(dmf[:,4]-dmf[:,1])/max(dmf[:,4]) * 100, "g--", linewidth=2)
# ax[0,0].set_xlim(11000,11500)
# ax[1,0].set_xlim(11000,11500)
# ax[2,0].set_xlim(11000,11500)
# ax[1,0].set_xlim(400,1300)
# ax[2,0].set_xlim(400,1300)
ax[1,0].set_ylabel('Vs Yspec(%)', fontsize=BIGGER_SIZE)
ax[2,0].set_ylabel('Vs Mineos(%)', fontsize=BIGGER_SIZE)
# ax[1,0].legend(['Mine vs Yspec', 'Mine vs Mineos'], loc='upper right', fontsize=BIGGER_SIZE)

######################################################

ax[0,1].plot(dmf[:,0], dmf[:,8], "b", linewidth=2)
ax[0,1].plot(dmf[:,0], dmf[:,5], "r--", linewidth=2)
ax[0,1].plot(dmf[:,0], dmf[:,2], "g--", linewidth=2)
ax[0,1].set_xlim(400,1300)


ax[1,1].plot(dmf[:,0],(dmf[:,5]-dmf[:,8])/max(dmf[:,5]) * 100, "k", linewidth=2)
ax[2,1].plot(dmf[:,0],(dmf[:,5]-dmf[:,2])/max(dmf[:,5]) * 100, "g--", linewidth=2)
ax[1,1].set_xlim(400,1300)
ax[2,1].set_xlim(400,1300)
# ax[1,1].plot(dmf[:,0]/3600,(dmf[:,5]-dmf[:,2])/max(dmf[:,5]) * 100, "b--", linewidth=2)

######################################################
# ax[0,2].plot(d[:,0]/3600,d[:,3], "g--", linewidth=3)
# ax[0,2].plot(d1[:,0]/3600, d1[:,3]  , "b", linewidth=2)
# ax[0,2].plot(d2[:,0]/3600, d2[:,3], "r--", linewidth=2)
# ax[0,2].plot(dmf[:,0]/3600, dmf[:,3], "g--", linewidth=2)

ax[0,2].plot(dmf[:,0], dmf[:,9], "b", linewidth=2)
ax[0,2].plot(dmf[:,0], dmf[:,6] , "r--", linewidth=2)
ax[0,2].plot(dmf[:,0], dmf[:,3], "g--", linewidth=2)
ax[0,2].set_xlim(400,1300)


ax[1,2].plot(dmf[:,0],(dmf[:,6]-dmf[:,9])/max(dmf[:,6]) * 100, "k", linewidth=2)
ax[2,2].plot(dmf[:,0],(dmf[:,6]-dmf[:,3])/max(dmf[:,6]) * 100, "g--", linewidth=2)
ax[1,2].set_xlim(400,1300)
ax[2,2].set_xlim(400,1300)
# ax[1,2].plot(dmf[:,0]/3600,(dmf[:,6]-dmf[:,3])/max(dmf[:,6]) * 100, "b--", linewidth=2)


######################################################


ax[0,0].set_title('Z', fontsize=BIGGER_SIZE)
ax[0,1].set_title('North', fontsize=BIGGER_SIZE)
ax[0,2].set_title('East', fontsize=BIGGER_SIZE)

ax[1,0].set_xlabel('Time (hr)', fontsize=BIGGER_SIZE)
ax[1,1].set_xlabel('Time (hr)', fontsize=BIGGER_SIZE)
ax[1,2].set_xlabel('Time (hr)', fontsize=BIGGER_SIZE) 
plt.show()