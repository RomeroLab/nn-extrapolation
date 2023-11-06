# Neural network extrapolation to distant regions of the protein fitness landscape


## Installation
The code in this module is built on the [nn4dms](github.com/gitter-lab/nn4dms) machine learning models and module. This module has the same requirements as nn4dms, with a few additional dependencies. Use [conda](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#create-env-file-manually) to set up an environment from [env.yml](env.yml)

```
conda env create -f env.yml
conda activate gb1_inf
```

Installation takes ~5 minutes.

## Software Requirements
This module has been tested on the following systems:
- macOS: Mojave (10.14.6)
- Linux: Green Obsidian (8.8)

## Data Analysis
### Extrapolating Learned Protein Fitness Lanscapes

```
python 01_extrapolation_predictions.py
python 01_extrapolation_trajectories.py
```


