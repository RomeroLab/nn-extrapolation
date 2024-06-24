#!/usr/bin/env python
# coding: utf-8


import os
from os.path import abspath, join, isdir
import sys

# for relative paths in nn4dms code to work properly, we need to set the current working
# directory to the root of the project
# we also need to add the code folder to the system path for imports to work properly

print('Setting working directory to nn4dms root.')
os.chdir('nn4dms_nn-extrapolate')
module_path = abspath("code")
if module_path not in sys.path:
    sys.path.append(module_path)

# add relative path to write directory (nn-extrapolate)
nnextrap_root_relpath = '..'
pretrained_dir = "nn-extrapolation-models/pretrained_models"

import design_tools as tools
import pickle
import random
import numpy as np
import sys
import yaml
import importlib


AAs = 'ACDEFGHIKLMNPQRSTVWY'

def load_config(config_file):
    with open(join(nnextrap_root_relpath, config_file), 'r') as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

def run_simulated_annealing(config):
    AA_options = [tuple([AA for AA in AAs]) for i in range(len(config['WT']))]
    AA_options.pop(0)
    AA_options.insert(0, ['M'])

    seq2fitness_tools = importlib.__import__(config['seq2fitness_tools'])
    seq2fitness_handler = seq2fitness_tools.seq2fitness_handler()
    print('setting up optimizer...')
    sa_optimizer = tools.SA_optimizer(seq2fitness_handler.seq2fitness, config['WT'], AA_options,
            config['num_mut'], mut_rate=config['mut_rate'], nsteps=config['nsteps'],
            cool_sched=config['cool_sched'])
    print('running optimization...')
    best_mut, fitness = sa_optimizer.optimize(seed=config['seed'])
    with open(join(nnextrap_root_relpath, config['export_best_seqs']), 'wb') as f:
        pickle.dump([best_mut, fitness], f)

    if config['save_plot_trajectory']:
        sa_optimizer.plot_trajectory(savefig_name=join(nnextrap_root_relpath, config['file_plot_trajectory']))


if __name__ == '__main__':
    print('running 02_run_sa.py', sys.argv[1])
    config = load_config(sys.argv[1])
    run_simulated_annealing(config)

