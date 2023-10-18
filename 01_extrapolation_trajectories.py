#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from os.path import isfile, join, exists
from os import listdir
from scipy.stats import spearmanr


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


CHARS = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
         "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
WT = "MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE"

models = ['lr', 'fcn', 'gcn', 'cnn']

for direction in ['all', 'all_down', 'wt']:
    func_all_list = []
    label_list = []

    for model in models:
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

        # calculate wt fitness for each model
        if direction == 'wt':
            seqs = [WT]
            func_all = []
            encoded_variants = enc.encode(encoding="one_hot,aa_index", char_seqs=seqs, wt_aa=[aa for aa in WT])

            functions_all = []
            for sess in model_sesses: # TODO: change back to all cnns
                functions_all.append(inf.run_inference(encoded_data=encoded_variants, sess=sess))
            functions = np.median(functions_all, axis=0)
            seqs_mut_df = pd.DataFrame(data=list(zip(seqs, np.array(functions_all).T, [functions])), 
                                       columns=['seq', 'func_all', 'func'])
            func_all_list.append([np.array(functions_all)])
            label_list.append(model+"_func")
            
        # make upward or downward trajectory for each model
        else:
            curr_muts = []
            func_all = []
            med_fits = []
            num_muts = 55

            for i in range(55): 
                print("starting mutation ", i)
                curr_poss = [int(mut[1:-1]) for mut in curr_muts]

                # look at all possible neighbors
                possible_muts = []
                for pos in range(1, len(WT)):
                    if pos not in curr_poss: # don't mutate already mutated positions
                        for aa in CHARS:
                            if WT[pos] != aa: # not mutating to self
                                mut = WT[pos]+str(pos)+aa
                                possible_muts.append(mut)

                # calculate fitness
                join_muts = [curr_muts+[mut] for mut in possible_muts]
                seqs = [dt.mut2seq(WT, muts) for muts in join_muts]
                encoded_variants = enc.encode(encoding="one_hot,aa_index", char_seqs=seqs, wt_aa=[aa for aa in WT])

                functions_all = []
                for sess in model_sesses: 
                    functions_all.append(inf.run_inference(encoded_data=encoded_variants, sess=sess))
                functions = np.median(functions_all, axis=0)
                # get next mutant in trajectory as min/max of median model prediction
                seqs_mut_df = pd.DataFrame(data=list(zip(possible_muts, join_muts, seqs, np.array(functions_all).T, functions)), 
                                           columns=['added_mut', 'join_mut', 'seq', 'func_all', 'func'])
                if direction == 'all'
                    seqs_mut_df.sort_values('func', ascending=False, inplace=True)
                else:
                    seqs_mut_df.sort_values('func', ascending=True, inplace=True)
                curr_muts.append(seqs_mut_df.iloc[0]['added_mut'])
                func_all.append(seqs_mut_df.iloc[0]['func_all'])
                med_fits.append(seqs_mut_df.iloc[0]['func'])
                num_muts -= 1
                if num_muts <= 0:
                    break

            # save steps
            func_all_list.append(curr_muts)
            func_all_list.append(func_all)
            label_list.append(model+"_mut")
            label_list.append(model+"_func")


    func_df = pd.DataFrame(data=list(zip(*func_all_list)), columns=label_list)
    func_df.to_csv('gen_data/mut_func_'+direction+'.csv')
