"""Refine time quadrature of the identical saved, filtered contour spectrum."""
from pathlib import Path
import argparse
import json
import numpy as np

p = argparse.ArgumentParser()
p.add_argument('directory', type=Path)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
s = np.genfromtxt(a.directory / 'spectrum.csv', delimiter=',', names=True)
df = float(np.diff(s['frequency_hz'])[0])
f = s['frequency_hz']
b = np.where(f < .008, 1., np.where(f < .01, .5 * (1 + np.cos(np.pi * (f - .008) / .002)), 0.))
values = (s['real'] + 1j * s['imag']) * b
rows = []
for dt in [5., 2.5, 1.25, .625]:
    n = round(1 / (df * dt))
    z = np.zeros(n // 2 + 1, dtype=complex)
    z[s['bin'].astype(int)] = values
    t = np.arange(n) * dt
    u = np.fft.irfft(z, n=n) * n * df * np.exp((5 / 3600) * t)
    window = np.where((t >= 1000) & (t <= 1200), .5 * (1 - np.cos(2 * np.pi * (t - 1000) / 200)), 0.)
    rows.append({'dt_s': dt, 'n': n, 'energy_m2_s': float(.5 * dt * np.sum(window * u**2))})
for row in rows:
    row['relative_to_finest_energy'] = abs(row['energy_m2_s'] / rows[-1]['energy_m2_s'] - 1)
result = {'spectrum_directory': str(a.directory), 'df_hz': df,
          'policy': 'Fixed saved contour coefficients, epsilon, filter and FFT period; zero-pad high frequency bins to refine only time sampling and Hann quadrature. No additional Earth solve.',
          'rows': rows, 'energy_gate_1e_3_passed': rows[0]['relative_to_finest_energy'] <= 1e-3}
a.output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
