#!/usr/bin/env python3
"""Descriptive degree comparison of saved sections; no new wave solves."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.signal import resample
from check_phase_b import save

p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('geometry',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
sections={};metadata={}
for degree in [384,1536]:
    folder=a.root/f'section-l{degree}/data'
    m=json.loads((folder/'section.json').read_text());metadata[degree]=m
    angles=np.genfromtxt(folder/'section_angles.csv',delimiter=',',names=True)
    radial=np.genfromtxt(folder/'section_radial_sides.csv',delimiter=',',names=True)
    values=np.memmap(folder/'section.bin',dtype='f8',mode='r').reshape(179,m['angles'],6)
    source=int(np.argmin(np.abs(angles['along_deg'])))
    assert angles['along_deg'][source]==0
    sections[degree]=np.roll(values,-source,axis=1)
fine=sections[1536];n=fine.shape[1];along=360*np.arange(n)/n
receiver=json.loads(a.geometry.read_text())['receiver_along_deg']
distance=lambda point:np.abs((along-point+180)%360-180)
mask=(distance(0)>=5)&(distance(receiver)>=5)
rows=[]
for c,name in enumerate(['A','C','F','L','N','kappa']):
    coarse=resample(sections[384][:,:,c],n,axis=1)
    x=coarse[:,mask];y=fine[:,:,c][:,mask]
    norm=float(np.linalg.norm(y));scale=float(np.max(np.abs(y)))
    error=x-y
    rows.append({'channel':name,'relative_l2_vs_fine':float(np.linalg.norm(error)/norm) if norm else 0,
                 'max_difference_scaled_by_fine_max':float(np.max(np.abs(error))/scale) if scale else 0})
energies={degree:json.loads((a.root/f'validate-l{degree}/checks.json').read_text())['energy_m2_s'] for degree in [384,1536]}
save(a.output,{'status':'DESCRIPTIVE_COMPARISON_NOT_CONVERGENCE_CERTIFICATE',
               'degrees':[384,1536],'channels':rows,'radial_sides':179,
               'fine_section_angles':n,'retained_angles_outside_5degree_masks':int(mask.sum()),
               'alignment':'periodic Fourier resampling of the complete bandlimited L384 section onto L1536 section angles; source at0deg',
               'energy_m2_s':energies,'energy_difference_relative_to_fine':abs(energies[1536]-energies[384])/abs(energies[1536]),
               'limits':['Two degrees only; no established convergence.','No radial interpolation; both interface sides remain distinct.',
                         'Compares raw energy derivatives, with degree-dependent waveforms/adjoint weights.','No YSpec, full volume or new radial solves.']})
