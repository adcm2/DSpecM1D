#%%
import numpy as np
import matplotlib.pyplot as plt

# --- 1. Load data_3 ---
# Use a try-except block for better error handling.
try:
    path_3 = "bench_step_error_3.out"
    data_3 = np.loadtxt(path_3, delimiter=";")
    path_10 = "bench_step_error_10.out"
    data_10 = np.loadtxt(path_10, delimiter=";")
    path_30 = "bench_step_error_30.out"
    data_30 = np.loadtxt(path_30, delimiter=";")
    path_50 = "bench_step_error_50.out"
    data_50 = np.loadtxt(path_50, delimiter=";")

except FileNotFoundError:
    print(f"Error: data_3 file not found at '{path_3}'")
    exit()

# --- 2. Extract data_3 Columns ---
# Give meaningful names to the columns being plotted.
step_sizes = data_3[:, 0]
# Assuming the columns are: time-domain error (z, n), freq-domain error (z, n)
t_error_z = data_3[:, 1]
t_error_n = data_3[:, 2]
t_error_e = data_3[:, 3]
w_error_z = data_3[:, 4]
w_error_n = data_3[:, 5]
w_error_e = data_3[:, 6]

t_error_1 = (t_error_z + t_error_n + t_error_e) / 3.0
w_error_1 = (w_error_z + w_error_n + w_error_e) / 3.0

# do the same for data_10 if needed
step_sizes_10 = data_10[:, 0]
t_error_10 = (data_10[:, 1] + data_10[:, 2] + data_10[:, 3]) / 3.0
w_error_10 = (data_10[:, 4] + data_10[:, 5] + data_10[:, 6]) / 3.0

# do the same for data_30 if needed
step_sizes_30 = data_30[:, 0]
t_error_30 = (data_30[:, 1] + data_30[:, 2] + data_30[:, 3]) / 3.0
w_error_30 = (data_30[:, 4] + data_30[:, 5] + data_30[:, 6]) / 3.0
# do the same for data_50 if needed
step_sizes_50 = data_50[:, 0]
t_error_50 = (data_50[:, 1] + data_50[:, 2] + data_50[:, 3]) / 3.0
w_error_50 = (data_50[:, 4] + data_50[:, 5] + data_50[:, 6]) / 3.0

# --- 3. Setup Figure ---
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
plt.style.use('seaborn-v0_8-whitegrid')

# --- 4. Define Plotting Parameters ---
lwidth = 2.0
BIGGER_SIZE = 16
marker_style = 'o'
ms_t = 'o'
ms_w = 's'

# --- 5. Plot data_3 ---
# Plot the error for each component on a log-log scale.
# ax.plot(step_sizes, t_error_z, "b-", marker=marker_style, linewidth=lwidth, label='Time-domain (Z)')
# ax.plot(step_sizes, t_error_n, "r-", marker=marker_style, linewidth=lwidth, label='Time-domain (N)')
# ax.plot(step_sizes, t_error_e, "g-", marker=marker_style, linewidth=lwidth, label='Time-domain (E)')
# ax.plot(step_sizes, w_error_z, "b--", marker=marker_style, linewidth=lwidth, label='Freq-domain (Z)')
# ax.plot(step_sizes, w_error_n, "r--", marker=marker_style, linewidth=lwidth, label='Freq-domain (N)')
# ax.plot(step_sizes, w_error_e, "g--", marker=marker_style, linewidth=lwidth, label='Freq-domain (E)')

# ax.plot(step_sizes, t_error_1, "b-", marker=ms_t, linewidth=lwidth , label='Time-domain (Avg)')
# ax.plot(step_sizes, w_error_1, "b--", marker=ms_w, linewidth=lwidth, label='Freq-domain (Avg)')
# ax.plot(step_sizes_10, t_error_10, "r-", marker=ms_t, linewidth=lwidth, label='Time-domain (Avg) 10mHz')
# ax.plot(step_sizes_10, w_error_10, "r--", marker=ms_w, linewidth=lwidth, label='Freq-domain (Avg) 10mHz')

# ratio of step size to min wavelength:
ratio_3 = 0.63/(3 * step_sizes)
ratio_10 = 0.63/(10 * step_sizes_10)
ratio_30 = 0.63/(30 * step_sizes_30)
ratio_50 = 0.63/(50 * step_sizes_50)
ax.plot(ratio_3, t_error_1, "b-",linewidth=lwidth , label='3mHz')
# ax.plot(ratio_3, w_error_1, "b--",  linewidth=lwidth, label='3mHz, freq-domain')
ax.plot(ratio_10, t_error_10, "r-",linewidth=lwidth, label='10mHz')
# ax.plot(ratio_10, w_error_10, "r--", linewidth=lwidth, label='10mHz, freq-domain')
ax.plot(ratio_30, t_error_30, "g-",linewidth=lwidth, label='30mHz')
ax.plot(ratio_50, t_error_50, "k-",linewidth=lwidth, label='50mHz')
# ax.plot(ratio_50, w_error_50, "g--", linewidth=lwidth, label='50mHz, freq-domain')

# --- 6. Final Figure Adjustments ---
# Set log-log scale
ax.set_xscale('log')
ax.set_yscale('log')

# Set labels and title
ax.set_xlabel('# elements per wavelength', fontsize=BIGGER_SIZE)
ax.set_ylabel('Average relative error (%)', fontsize=BIGGER_SIZE)
# ax.set_title('Convergence with Timestep Size', fontsize=BIGGER_SIZE + 2)

# Customize ticks and legend
ax.tick_params(axis='both', which='major', labelsize=BIGGER_SIZE - 2)
ax.legend(fontsize=BIGGER_SIZE - 2)

# Improve layout
plt.tight_layout()
plt.show()