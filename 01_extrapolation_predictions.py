#!/usr/bin/env python
# coding: utf-8

# ## Calculate predictions for Wu combinatorial Library
# Ran this section on the group server

import numpy as np
import pandas as pd
from os.path import isfile, join, exists
from os import listdir

from os.path import abspath
import sys
module_path = abspath("nn4dms/code")
if module_path not in sys.path:
    sys.path.append(module_path)

import constants
import utils
import encode as enc
import inference as inf
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
df = pd.read_csv('data/elife-16965-supp1-v4.csv')


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
encoded_variants = enc.encode(encoding="one_hot,aa_index", char_seqs=df.sequence.tolist(), wt_aa=[aa for aa in WT])

models = ['lr', 'fcn', 'gcn', 'cnn']

ind_model_path = 'pretrained_models/other_models/gb1_'
for model in models:
    # model_sess = inf.restore_sess(path + model)
    with inf.restore_sess(ind_model_path + model) as model_sess:
        df[model+'_pred'] = inf.run_inference(encoded_data=encoded_variants, sess=model_sess)

    num_models = 100
    found_models = []
    model_paths = []
    for i in range(num_models):
        path = 'pretrained_models/'+model+'s/model_'+str(i)
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

    # run inferences
    model_pred_all = []
    for sess in model_sesses:
        model_pred_all.append(inf.run_inference(encoded_data=encoded_variants, sess=sess))
    
    df[model+'_pred_all'] = list(zip(*model_pred_all))

# save to csv
df.to_csv('gen_data/pred_extrapolation_wu.csv')
