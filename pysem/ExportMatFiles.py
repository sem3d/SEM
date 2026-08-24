# -*- coding: utf-8 -*-
"""
Le os traces H5 de uma simulacao SEM (capteurs), salva a energia global como
no GetCapteurs.py, e exporta cada variavel presente no H5 ('Variables') no seu
proprio arquivo .mat (nome em ingles), sem reamostrar em nenhuma grade
particular. Nao inclui o calculo dos corretores de homogeneizacao (isso fica
em GetCapteurs.py) -- so a leitura/exportacao generica dos traces.

Funciona tanto para SEM2D quanto SEM3D: o layout de colunas vem direto do
dataset 'Variables' de cada arquivo H5, entao o numero de componentes por
campo (2D: x,z; 3D: x,y,z) nao precisa ser conhecido de antemao.
"""

import itertools
import os
import re

import h5py
import numpy as np
import hdf5storage

# Diretorio com os capteurs*.h5, arquivo de estacoes e destino dos .mat -- edite aqui.
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


# 'Variables' e diretamente o cabecalho das colunas de cada dataset de estacao: uma
# entrada por componente, no formato "<nome><espacos><indice da componente>" (ex.:
# "Displ      1", "Displ      2", "Displ      3"). O layout dos campos (nome, largura)
# sai direto dai -- sem precisar manter uma lista hardcoded sincronizada a mao toda vez
# que o simulador muda as variaveis de saida.
FIELD_LABEL_RE = re.compile(r'^(.*\S)\s+\d+$')


def parse_field_layout(labels):
    names = [FIELD_LABEL_RE.match(lab).group(1).replace(' ', '') if FIELD_LABEL_RE.match(lab) else lab.replace(' ', '')
             for lab in labels]
    return [(name, len(list(group))) for name, group in itertools.groupby(names)]


# Traducao em ingles dos nomes de campo, usada para nomear os .mat de saida.
# Um campo ausente deste dicionario mantem o nome original (ver mais abaixo).
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

# Descobre o layout ('Variables') e o numero de passos de tempo (apos subamostragem)
# sem ler todos os dados -- so olha metadados do primeiro arquivo que tiver o necessario.
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
    raise RuntimeError("'Variables' nao encontrado, ou nenhum dataset de estacao nos arquivos H5")

# coun comeca em 0: a coluna 0 dos dados JA e a primeira variavel declarada em
# 'Variables' (tipicamente 'Time'), nao uma coluna implicita antes dela.
offsets = {}
coun = 0
for field_name, width in parse_field_layout(labels):
    offsets[field_name] = (coun, width)
    coun += width

# Le direto em arrays indexados pelo numero da estacao (0-based, o N em "UU_N"),
# sem passar por uma lista intermediaria de dicts nem por um sort no final.
N = line_count
Time = np.empty((N, num_time))
fields = {name: np.empty((N, num_time, width) if width > 1 else (N, num_time))
          for name, (offset, width) in offsets.items()}
filled = np.zeros(N, dtype=bool)
DataE = []

capcount = 1
for fname in files:
    with h5py.File(os.path.join(traces_dir, fname), 'r') as f:
        print("Lendo", fname)
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
            # .copy() mantido de proposito: sem ele, algumas estacoes voltaram com
            # tamanho inconsistente em execucoes anteriores.
            aux = f[dataset_name][::subs, :].copy()
            Time[row] = aux[:, 0]
            for field_name, (offset, width) in offsets.items():
                fields[field_name][row] = aux[:, offset] if width == 1 else aux[:, offset:offset + width]
            filled[row] = True

missing = np.flatnonzero(~filled)
if missing.size:
    raise RuntimeError(f"Faltam {missing.size} estacoes (ids 0-based, ex.: {missing[:10].tolist()})")

ttime = Time[0]

os.makedirs(out_dir, exist_ok=True)

# hdf5storage.write() nao sobrescreve de forma limpa um .mat corrompido/incompleto de uma
# rodada anterior (falha com "bad object header version number") -- apaga antes de escrever.
output_names = ['VarEnergia.mat'] + [f'{FIELD_TO_ENGLISH.get(name, name)}.mat' for name in fields]
for output_name in output_names:
    try:
        os.remove(os.path.join(out_dir, output_name))
    except FileNotFoundError:
        pass

if DataE:
    print("Salvando a energia")
    hdf5storage.write({'E': DataE}, '.', os.path.join(out_dir, 'VarEnergia.mat'), matlab_compatible=True)

# Um arquivo .mat por variavel presente no H5, nomeado em ingles
for field_name, data in fields.items():
    english_name = FIELD_TO_ENGLISH.get(field_name, field_name)
    print(f"Salvando {english_name}.mat")
    hdf5storage.write(
        {'Time': ttime, english_name: data},
        '.', os.path.join(out_dir, f'{english_name}.mat'),
        matlab_compatible=True,
    )

print("Exportacao concluida!")
