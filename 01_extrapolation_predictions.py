#!/usr/bin/env python
# coding: utf-8

# ## Calculate predictions for Wu combinatorial Library
# Ran this section on the group server

import numpy as np
import pandas as pd
from os.path import isfile, join, exists
from os import listdir, chdir
from tqdm import tqdm

from os.path import abspath
import sys


# for relative paths in nn4dms code to work properly, we need to set the current working
# directory to the root of the project
# we also need to add the code folder to the system path for imports to work properly
print('Setting working directory to nn4dms root.')
chdir('nn4dms_nn-extrapolate')
module_path = abspath("code")
if module_path not in sys.path:
    sys.path.append(module_path)

# add relative path to write directory (nn-extrapolation)
nnextrap_root_relpath = ".."
pretrained_dir = "nn-extrapolation-models/pretrained_models"

import constants
import utils
import encode as enc
import inference as inf
import inference_lr as inf_lr
import design_tools as dt

''' 
Calculate enrichment for a list of variants
Args:
    wt_unsel: WT count in unsorted population
    wt_sel: WT count in sorted population
    vars_unsel: numpy array of variant counts in unsorted population
    vars_sel: numpy array of variant counts in sorted population
Returns:
    numpy array containing enrichment scores for all variants
'''
def calc_enrich(wt_unsel, wt_sel, vars_unsel, vars_sel):
    assert wt_unsel > 0
    assert wt_sel > 0
    return np.log((vars_sel+0.5)/(wt_sel+0.5)) - np.log((vars_unsel+0.5)/(wt_unsel+0.5))

CHARS = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
         "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
WT = "MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE"

AAs = CHARS

positions = [38, 39, 40, 53]
WT_res = "VDGV" # TODO: Check that this is WT!

# load data from Wu et al.
df = pd.read_csv(join(nnextrap_root_relpath, 'data/elife-16965-supp1-v4.csv'))


# get full amino acid sequnce
sequences = []
for var in df.Variants:
    sequence = WT
    for wt_res, pos, res in zip(WT_res, positions, var):
        assert(wt_res == sequence[pos])
        sequence = sequence[:pos] + res + sequence[pos+1:]
        
    sequences.append(sequence)
df['sequence'] = sequences

# calc enrich2 fitnesses
wt_data = df.loc[df.Variants == WT_res]
wt_unsel = wt_data['Count input'].values[0]
wt_sel = wt_data['Count selected'].values[0]

df['enrich2_fit'] = calc_enrich(wt_unsel, wt_sel, df['Count input'].to_numpy(), df['Count selected'].to_numpy())


# encode variants
print('Encoding variants...')
encoded_variants = enc.encode(encoding="one_hot,aa_index", char_seqs=df.sequence.tolist(), wt_aa=[aa for aa in WT])

models = ['lr', 'fcn', 'gcn', 'cnn']

ind_model_path = join(nnextrap_root_relpath, pretrained_dir, 'other_models/gb1_')
print('Calculating fitnesses for LR, GCN, GCN, and CNN models...')
for model in tqdm(models, total=len(models), ncols=100, desc="Model"):
    # get fitnesses from individual models used in paper
    # lr requires separate inference to remove ph parameter
    with inf.restore_sess(ind_model_path + model) as model_sess:
        if model == 'lr':
            df[model+'_pred'] = inf_lr.run_inference_lr(encoded_data=encoded_variants, sess=model_sess)
        # use inf import for all other models
        else:
            df[model+'_pred'] = inf.run_inference(encoded_data=encoded_variants, sess=model_sess)

    num_models = 100
    found_models = []
    model_paths = []
    for i in tqdm(range(num_models), total=num_models, ncols=100, leave=False):
        path = join(nnextrap_root_relpath, pretrained_dir, model+'s/model_'+str(i))
        if (exists(path)):
            for file_name in listdir(path):
                if '.pb' in file_name:
                    model_name = file_name
                    model_paths.append(path+'/'+model_name)
                    found_models.append(i)
    if len(model_paths) != num_models:
        print('Could not find all models, missing models: ',
                ','.join([str(i) for i in range(100) if i not in found_models]))
    model_sesses = []
    for model_path in model_paths:
        model_sesses.append(inf.restore_sess_from_pb(model_path))

    # run inferences for additional models
    model_pred_all = []
    for sess in model_sesses:
        if model == 'lr':
            model_pred_all.append(inf_lr.run_inference_lr(encoded_data=encoded_variants, sess=sess))
        else:
            model_pred_all.append(inf.run_inference(encoded_data=encoded_variants, sess=sess))
    
    df[model+'_pred_all'] = list(zip(*model_pred_all))


# further process predictions to save only essential data used for processing. Must do this because full
# output is ~1.3 GB, which is too big for GitHub.

# get fitness for EnsC and EnsM by taking 5th lowest value and median value from cnn_pred_all
df['cnn_pred_all_sort'] = [pd.Series(ele).sort_values() for ele in df['cnn_pred_all']]
df['ensc_pred'] = [ele.iloc[5] for ele in df['cnn_pred_all_sort']]
df['ensm_pred'] = [ele.median() for ele in df['cnn_pred_all_sort']]

# columns to save
save_columns = [
    'Variants',
    'HD',
    'enrich2_fit',
    'lr_pred',
    'fcn_pred',
    'cnn_pred',
    'ensc_pred',
    'ensm_pred',
    'gcn_pred',
]

# save processed data to csv
df[save_columns].to_csv(join(nnextrap_root_relpath, 'gen_data/pred_extrapolation_wu.csv'))
# save full data to csv
df.to_csv(join(nnextrap_root_relpath, 'gen_data/raw_pred_extrapolation_wu.csv'))
