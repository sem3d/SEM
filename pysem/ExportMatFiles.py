# -*- coding: utf-8 -*-
"""
Reads H5 traces from a SEM simulation (capteurs/receivers), saves the global energy
similarly to GetCapteurs.py, and exports each variable present in H5 ('Variables') into
its own .mat file (English name), without resampling onto any particular grid.
Does not include homogenization corrector calculations (handled in GetCapteurs.py)
-- only generic trace reading/export.

Works for both SEM2D and SEM3D: the column layout comes directly from the
'Variables' dataset in each H5 file, so the number of components per
field (2D: x,z; 3D: x,y,z) does not need to be known in advance.
"""

import itertools
import os
import re

import h5py
import numpy as np
import hdf5storage

# Directory with capteurs*.h5 files, stations file, and .mat output destination -- edit here.
traces_dir = './traces/'
stations_path = os.path.join(traces_dir, os.pardir, 'stations.txt')
out_dir = traces_dir
subs = 10


def parse_index(dataset_name):
    parts = dataset_name.split('_')
    try:
        return int(parts[1])
    except (IndexError, ValueError):
        return None


# 'Variables' is directly the column header of each station dataset: one
# entry per component, in the format "<name><spaces><component index>" (e.g.,
# "Displ      1", "Displ      2", "Displ      3"). The field layout (name, width)
# comes directly from it -- without needing to manually maintain a hardcoded
# list whenever the simulator changes output variables.
FIELD_LABEL_RE = re.compile(r'^(.*\S)\s+\d+$')


def parse_field_layout(labels):
    names = [FIELD_LABEL_RE.match(lab).group(1).replace(' ', '') if FIELD_LABEL_RE.match(lab) else lab.replace(' ', '')
             for lab in labels]
    return [(name, len(list(group))) for name, group in itertools.groupby(names)]


# English translation of field names, used to name the output .mat files.
# Any field missing from this dictionary keeps its original name (see below).
FIELD_TO_ENGLISH = {
    'EnergyP': 'PotentialEnergy',
    'EnergyK': 'KineticEnergy',
    'EpsVol': 'VolumetricStrain',
    'Displ': 'Displacement',
    'Veloc': 'Velocity',
    'Accel': 'Acceleration',
    'Pressure': 'Pressure',
    'EpsDev': 'DeviatoricStrain',
    'EpsDevPl': 'PlasticDeviatoricStrain',
    'StressDev': 'DeviatoricStress',
    'EnergyD': 'DecomposedEnergy',
    'DUDX': 'DisplacementGradient',
    'GradLambda': 'LambdaGradient',
    'GradMu': 'MuGradient',
    'Rotat': 'Rotation',
}

with open(stations_path, 'r', encoding='utf-8') as f:
    line_count = sum(1 for _ in f)

files = sorted(f for f in os.listdir(traces_dir) if f.startswith('capteurs') and f.endswith('.h5'))

# Discover the layout ('Variables') and the number of time steps (after subsampling)
# without reading all data -- only inspects metadata of the first file containing them.
labels = None
num_time = None
for fname in files:
    with h5py.File(os.path.join(traces_dir, fname), 'r') as f:
        if labels is None and 'Variables' in f:
            labels = [lab.decode('utf-8').strip() for lab in f['Variables'][:]]
        if num_time is None:
            for dataset_name in f.keys():
                if dataset_name in ('Variables', 'Energy_Variables', 'Energy') or dataset_name.endswith('_pos'):
                    continue
                if parse_index(dataset_name) is not None:
                    num_time = len(range(0, f[dataset_name].shape[0], subs))
                    break
    if labels is not None and num_time is not None:
        break

if labels is None or num_time is None:
    raise RuntimeError("'Variables' not found, or no station dataset in H5 files")

# coun starts at 0: column 0 of the data is ALREADY the first variable declared in
# 'Variables' (typically 'Time'), not an implicit column before it.
offsets = {}
coun = 0
for field_name, width in parse_field_layout(labels):
    offsets[field_name] = (coun, width)
    coun += width

# Read directly into arrays indexed by station number (0-based, the N in "UU_N"),
# without going through an intermediate list of dicts or sorting at the end.
N = line_count
Time = np.empty((N, num_time))
fields = {name: np.empty((N, num_time, width) if width > 1 else (N, num_time))
          for name, (offset, width) in offsets.items()}
filled = np.zeros(N, dtype=bool)
DataE = []

capcount = 1
for fname in files:
    with h5py.File(os.path.join(traces_dir, fname), 'r') as f:
        print("Reading", fname)
        for dataset_name in sorted(f.keys()):
            if dataset_name in ('Variables', 'Energy_Variables') or dataset_name.endswith('_pos'):
                continue
            if dataset_name == 'Energy':
                DataE.append(f[dataset_name][:])
                continue
            cap_index = parse_index(dataset_name)
            if cap_index is None:
                continue
            print(f'{capcount}/{line_count}')
            capcount += 1
            row = cap_index
            # .copy() kept intentionally: without it, some stations returned with
            # inconsistent sizes in previous runs.
            aux = f[dataset_name][::subs, :].copy()
            Time[row] = aux[:, 0]
            for field_name, (offset, width) in offsets.items():
                fields[field_name][row] = aux[:, offset] if width == 1 else aux[:, offset:offset + width]
            filled[row] = True

missing = np.flatnonzero(~filled)
if missing.size:
    raise RuntimeError(f"Missing {missing.size} stations (0-based IDs, e.g.: {missing[:10].tolist()})")

ttime = Time[0]

os.makedirs(out_dir, exist_ok=True)

# hdf5storage.write() does not cleanly overwrite a corrupted/incomplete .mat from a
# previous run (fails with "bad object header version number") -- remove before writing.
output_names = ['ModelEnergy.mat'] + [f'{FIELD_TO_ENGLISH.get(name, name)}.mat' for name in fields]
for output_name in output_names:
    try:
        os.remove(os.path.join(out_dir, output_name))
    except FileNotFoundError:
        pass

if DataE:
    print("Saving energy")
    hdf5storage.write({'E': DataE}, '.', os.path.join(out_dir, 'ModelEnergy.mat'), matlab_compatible=True)

# One .mat file per variable present in H5, named in English
for field_name, data in fields.items():
    english_name = FIELD_TO_ENGLISH.get(field_name, field_name)
    print(f"Saving {english_name}.mat")
    hdf5storage.write(
        {'Time': ttime, english_name: data},
        '.', os.path.join(out_dir, f'{english_name}.mat'),
        matlab_compatible=True,
    )

print("Export completed!")
