#%%
import numpy as np
import matplotlib.pyplot as plt

# --- 1. Load data_3 ---
# Use a try-except block for better error handling.
try:
    path_3 = "../outputs/ex5_NQ5_step_error_50.out"
    data_3 = np.loadtxt(path_3, delimiter=";")
    path_6 = "../outputs/ex5_NQ6_step_error_50.out"
    data_6 = np.loadtxt(path_6, delimiter=";")
    path_4 = "../outputs/ex5_NQ4_step_error_50.out"
    data_4 = np.loadtxt(path_4, delimiter=";")

except FileNotFoundError:
    print(f"Error: data_3 file not found at '{path_3}'")
    exit()

# --- 2. Extract data_3 Columns ---
# Give meaningful names to the columns being plotted.
step_sizes_NQ5 = data_3[:, 0]
# Assuming the columns are: time-domain error (z, n), freq-domain error (z, n)
t_error_1 = (data_3[:, 1] + data_3[:, 2] + data_3[:, 3]) / 3.0
w_error_1 = (data_3[:, 4] + data_3[:, 5] + data_3[:, 6]) / 3.0

# NQ6
step_sizes_NQ6 = data_6[:, 0]
t_error_2 = (data_6[:, 1] + data_6[:, 2] + data_6[:, 3]) / 3.0
w_error_2 = (data_6[:, 4] + data_6[:, 5] + data_6[:, 6]) / 3.0


# NQ4
step_sizes_NQ4 = data_4[:, 0]
t_error_3 = (data_4[:, 1] + data_4[:, 2] + data_4[:, 3]) / 3.0
w_error_3 = (data_4[:, 4] + data_4[:, 5] + data_4[:, 6]) / 3.0
# --- 3. Setup Figure ---
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
plt.style.use('seaborn-v0_8-whitegrid')

# --- 4. Define Plotting Parameters ---
lwidth = 2.0
BIGGER_SIZE = 24
marker_style = 'o'
ms_t = 'o'
ms_w = 's'

# Colorblind-friendly palette (Okabe-Ito)
COLORS = {
    'nq4': '#009E73',  # bluish green
    'nq5': '#0072B2',  # blue
    'nq6': '#D55E00'   # vermillion
}


# ratio of step size to min wavelength:
ratio_4 = 0.63/(50 * step_sizes_NQ4)
ratio_5 = 0.63/(50 * step_sizes_NQ5)
ratio_6 = 0.63/(50 * step_sizes_NQ6)
ax.plot(ratio_4, t_error_3, color=COLORS['nq4'], linestyle='-', linewidth=lwidth, label='NQ = 4')
ax.plot(ratio_5, t_error_1, color=COLORS['nq5'], linestyle='-', linewidth=lwidth, label='NQ = 5')
ax.plot(ratio_6, t_error_2, color=COLORS['nq6'], linestyle='-', linewidth=lwidth, label='NQ = 6')
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