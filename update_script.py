import re

file_path = '/home/adcm2/space/c++/DSpecM1D_Draft/plotting/functions/ex1.py'
with open(file_path, 'r') as f:
    code = f.read()

new_func = """
def get_modes_in_range(fmin, fmax):
    mode_files = [
        "/home/adcm2/space/c++/mineos/DEMO/MYEX/lf_prem_R",
        "/home/adcm2/space/c++/mineos/DEMO/MYEX/lf_prem_S",
        "/home/adcm2/space/c++/mineos/DEMO/MYEX/lf_prem_T"
    ]
    modes = []
    for filepath in mode_files:
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
                start_idx = 0
                for i, line in enumerate(lines):
                    if "mode" in line and "phs vel" in line and "w(mhz)" in line:
                        start_idx = i + 2
                        break
                if start_idx > 0:
                    for line in lines[start_idx:]:
                        if not line.strip():
                            continue
                        parts = line.split()
                        if len(parts) >= 5:
                            try:
                                n = parts[0]
                                type_mode = parts[1].upper()
                                l = parts[2]
                                freq = float(parts[4])
                                if fmin <= freq <= fmax:
                                    name = f'${{_{n}}}${type_mode}${{_{{{l}}}}}$'
                                    modes.append((name, freq))
                            except ValueError:
                                pass
        except FileNotFoundError:
            pass
    modes.sort(key=lambda x: x[1])
    return [m[0] for m in modes], [m[1] for m in modes]

"""

func_insert_pos = code.find('# --- 7. Inset Plot (Zoomed Region) ---')
code = code[:func_insert_pos] + new_func + code[func_insert_pos:]


replacement1 = """
names, xvalues = get_modes_in_range(zoom_fmin, zoom_fmax)
yvalues = []
for x in xvalues:
    idx = np.searchsorted(zoom_x, x)
    if idx < len(zoom_x):
        window = 3
        start_w = max(0, idx - window)
        end_w = min(len(zoom_x), idx + window + 1)
        y_val1 = np.max(zoom_y_yspec[start_w:end_w]) if len(zoom_y_yspec[start_w:end_w]) > 0 else 0
        y_val2 = np.max(zoom_y_trans[start_w:end_w]) if len(zoom_y_trans[start_w:end_w]) > 0 else 0
        y_val3 = np.max(zoom_y_specnm[start_w:end_w]) if len(zoom_y_specnm[start_w:end_w]) > 0 else 0
        yvalues.append(max(y_val1, y_val2, y_val3))
    else:
        yvalues.append(0.05)
"""

old_str1 = """names = [r'${}_0$S${}_2$', r'${}_2$S${}_1$', r'${}_0$S${}_3$']
xvalues = [0.3108, 0.4063, 0.4713]
yvalues = [0.0832, 0.0096, 0.0556]"""

code = code.replace(old_str1, replacement1)


replacement2 = """
names2, xvalues2 = get_modes_in_range(zoom2_fmin, zoom2_fmax)

yvalues2 = []
for x in xvalues2:
    idx = np.searchsorted(zoom2_x, x)
    if idx < len(zoom2_x):
        window = 3
        start_w = max(0, idx - window)
        end_w = min(len(zoom2_x), idx + window + 1)
        y_val1 = np.max(zoom2_y_yspec[start_w:end_w]) if len(zoom2_y_yspec[start_w:end_w]) > 0 else 0
        y_val2 = np.max(zoom2_y_trans[start_w:end_w]) if len(zoom2_y_trans[start_w:end_w]) > 0 else 0
        y_val3 = np.max(zoom2_y_specnm[start_w:end_w]) if len(zoom2_y_specnm[start_w:end_w]) > 0 else 0
        yvalues2.append(max(y_val1, y_val2, y_val3))
    else:
        yvalues2.append(0.05)
"""

old_str2 = """names2 = [r'${}_4$S${}_4$', r'${}_1$S${}_{11}$', r'${}_5$S${}_4$',r'${}_2$S${}_{10}$', r'${}_6$S${}_2$', r'${}_0$S${}_{16}$',r'${}_2$S${}_0$',r'${}_7$S${}_2$']
xvalues2 = [2.292460, 2.359480, 2.387861, 2.425544, 2.460181, 2.472586,2.513545 ,2.532304]
yvalues2 = [0.081, 0.57, 0.12, 0.35, 0.005, 0.4, 0.007, 0.026]"""

code = code.replace(old_str2, replacement2)


label_replacement2 = """
    if i % 2 == 0:  # Alternate label placement height for clear reading
        ax_ins2.text(xvalues2[i], yvalues2[i] + 0.01, name, fontsize=14, fontweight='bold', ha='center')
    else:  
        ax_ins2.text(xvalues2[i], yvalues2[i] + 0.02, name, fontsize=14, fontweight='bold', ha='center')
"""

old_label2 = """    if i in [0,1,2,3,5,6,7]:  # For the more prominent peaks, place label above the line
        ax_ins2.text(xvalues2[i], yvalues2[i] + 0.01, name, fontsize=14, fontweight='bold', ha='center')
    else:  # For smaller peaks, place label slightly above the line to avoid overlap
        ax_ins2.text(xvalues2[i]-0.01, yvalues2[i] + 0.015, name, fontsize=14, fontweight='bold', ha='center')
"""

code = code.replace(old_label2, label_replacement2)

with open(file_path, 'w') as f:
    f.write(code)

