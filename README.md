# Neural network extrapolation to distant regions of the protein fitness landscape


## Installation
The code in this module is built on the [nn4dms](github.com/gitter-lab/nn4dms) machine learning models and module. This module has the same requirements as nn4dms, with a few additional dependencies. Use [conda](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#create-env-file-manually) to set up an environment from [env.yml](env.yml)

The model code is in the linked repo [nn4dms](https://github.com/gitter-lab/nn4dms/tree/2b4fcfd6c6e90321f21fa3264f677d639f33ba83). This code is required to reproduce some of our inference results. Code and models to generate model predictions are included as the `nn4dms_nn-extrapolate` and`nn-extrapolate-models` submodules. Download times may be long if cloning submodules. We include download commands for nn-extrapolation with and without submodules.

```
# clone repo without submodules
gh repo clone RomeroLab/nn-extrapolation

# checkout submodules individually
git submodule update nn4dms_nn-extrapolate
git submodule update nn-extrapolation-models

# clone full repo including submodules
gh repo clone RomeroLab/nn-extrapolation -- --recurse-submodules
```

We provide two conda environments to reproduce our analysis. `gb1_inf` is used to run all Makefile commands. `gb1_notebook` is used to reproduce data anlysis in the Jupyter notebooks.
```
# download and install environments
conda env create -f env_inference.yml
conda env create -f env_notebook.yml

# activate conda environments
conda activate gb1_inf
conda activate gb1_notebook
```

Installation takes ~5 minutes.

## Software Requirements
This module has been tested on the following systems:
- macOS: Mojave (10.14.6)
- Linux: Green Obsidian (8.8)

## Data Analysis
### Extrapolating learned protein fitness lanscapes
Generate predictions for the Wu et al. 1-4 mutant GB1 fitness dataset.
``` bash
python 01_extrapolation_predictions.py
```

Generate extrapolation trajectories.
``` bash
python 01_extrapolation_trajectories.py
```

Generate plots for Fig 1 and Fig S1 in `01_extrapolation_analysis.ipynb`

### ML-guided protein design for deep exploration of the fitness landscape
Design sequences using `02_run_sa.py`. See example below. Each design can take minutes to hours, depending on the model; this can be accelerated by running on a GPU.
``` bash
python 02_run_sa.py data/config_example.txt
```

Generate plots for Fig 2 and Fig S2-3 in `02_designs_analysis.ipynb`

### Large-scale experimental characterization of ML designed GB1 variants
From the raw fastq files, preprocess the results and determine counts for each variant in the library. (fastq files will be downloadable from the SRA; save the fastq files in a directory `fastq_files`). This step can take hours to days.
``` bash
mkdir merged_reads
for d in fastq_files/ ; do
python 03_preprocessing.py fastq_files/${d} ${d:0:6} merged_reads/ designs.csv designs_counts.csv
```

Generate plots for Fig 3 and Fig S4-8 in `03_design_experimental_analysis.ipynb`

### ML designed GB1s show improved display and IgG binding
Generate plots for Fig 4 and Fig S9 in `03_design_experimental_analysis.ipynb`
