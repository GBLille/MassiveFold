![header](imgs/massivefold_logo.svg)

[![DOI](https://zenodo.org/badge/617429072.svg)](https://doi.org/10.5281/zenodo.13870060)

MassiveFold is a tool that allows to massively expand the sampling of structure predictions by improving the computation 
of [AlphaFold](https://github.com/google-deepmind/alphafold) based predictions.

It optimizes the parallelization of the structure inference by splitting the computation on CPU for alignments, running 
automatically batches of structure predictions on GPU, and gathering the results in one global output directory, with a 
global ranking and a variety of plots.

MassiveFold uses [AFmassive](https://github.com/GBLille/AFmassive), [ColabFold](https://github.com/sokrypton/ColabFold) 
or [AlphaFold3](https://github.com/google-deepmind/alphafold3) as inference engine; AFmassive is an updated version of Björn 
Wallner's [AFsample](https://github.com/bjornwallner/alphafoldv2.2.0/) that offers additional diversity parameters for 
massive sampling.

## MassiveFold: parallelize protein structure prediction
MassiveFold's design (see schematic below) is optimized for GPU cluster usage. It allows fast computation for massive 
sampling by automatically splitting a large run of numerous predictions into several jobs. Each of these individual 
jobs are computed on a single GPU node and their results are then gathered as a single output with each prediction 
ranked on a global level instead of the level of each individual job.  

This automatic splitting is also convenient for massive sampling on a single GPU server to manage jobs priorities.  

MassiveFold is only available with the **SLURM** workload manager (Simple Linux Utility for Resource Management) as 
it relies heavily on its features (job array, job dependency, etc...).

![header](imgs/massivefold_diagram.svg)

A run is composed of three steps:  
1. **alignment**: on CPU, sequence alignments is the first step (can be skipped if alignments are already computed)

2. **structure prediction**: on GPU, structure predictions follow the massive sampling principle. The total number 
of predictions is divided into smaller batches and each of them is distributed on a single GPU. These jobs wait for the 
alignment job to be over, if the alignments are not provided by the user.

3. **post_treatment**: on CPU, it finishes the job by gathering all batches outputs and produces plots with the 
[plots module](#massivefold_plots-output-representation) to visualize the run's performances. This job is executed only once 
all the structure predictions are over. 

## Installation

MassiveFold was developed to run massive sampling with [AFmassive](https://github.com/GBLille/AFmassive), 
[ColabFold](https://github.com/sokrypton/ColabFold) and [AlphaFold3](https://github.com/google-deepmind/alphafold3), and relies on them for its installation.

Follow the MassiveFold [installation guide](docs/installation.md).  
It details these steps of the MassiveFold installation:

1. [**Retrieve MassiveFold**](docs/installation.md#1-retrieve-massivefold)

2. [**Install MassiveFold**](docs/installation.md#2-install-massivefold)

3. [**Create header files**](docs/installation.md#3-create-header-files)

4. [**Set custom parameters**](docs/installation.md#4-set-custom-parameters)

## Uninstallation

To uninstall MassiveFold, remove the three conda environments (`massivefold`, `mf-afmassive-1.1.5`, `mf-colabfold-1.5.5` and 
`mf-alphafold-3.0.1`) and remove the `MassiveFold` folder you cloned. Make sure you copy all the files and folders you want 
to keep from the `output` and `log` directories somewhere else. 

## Usage

Usage section includes the most simple way to run MassiveFold with examples. For more detail on its functionning and other cases, see the [usage documentation](docs/usage.md).

Activate the conda environment, then launch MassiveFold.
```bash
conda activate massivefold
./run_massivefold.sh -s <SEQUENCE_PATH> -r <RUN_NAME> -p <NUMBER_OF_PREDICTIONS_PER_MODEL> -f <JSON_PARAMETERS_FILE> -t <TOOL> 
```

Example for AFmassive:
```bash
./run_massivefold.sh -s input/H1140.fasta -r afm_default -p 5 -f AFmassive_params.json
```
Example for ColabFold:
```bash
./run_massivefold.sh -s input/H1140.fasta -r cf_default -p 5 -f ColabFold_params.json
```
Example for AlphaFold3:
```bash
./run_massivefold.sh -s input/H1140.fasta -r af3_default -p 5 -f AlphaFold3_params.json
```
## massivefold_plots: output representation

Additionally to configurating the plots parameters inside MassiveFold JSON param file, the plotting module can also be used on an already produced MassiveFold (or AlphaFold2) output to evaluate visually its predictions.  

For more details on this usage, see [MassiveFold plots documentation](docs/plots.md).

## Troubleshooting

Some known issues were identified and can be prevented by following steps described in the [troublehshooting documentation](docs/troubleshooting.md).

## Citation

If you use MassiveFold in your work, please cite:

Raouraoua N. et al. **MassiveFold: unveiling AlphaFold’s hidden potential with optimized and parallelized massive 
sampling**. 2024. **_Nature Computational Science_**, DOI: 10.1038/s43588-024-00714-4,  
https://www.nature.com/articles/s43588-024-00714-4  

## Authors
Nessim Raouraoua (UGSF - UMR 8576, France)  
Claudio Mirabello (NBIS, Sweden)  
Thibaut Véry (IDRIS, France)  
Christophe Blanchet (IFB, France)  
Björn Wallner (Linköping University, Sweden)  
Marc F Lensink (UGSF - UMR8576, France)  
Guillaume Brysbaert (UGSF - UMR 8576, France)  

This work was carried out as part of Work Package 4 of the [MUDIS4LS project](https://www.france-bioinformatique.fr/actualites/mudis4ls-le-projet-despaces-numeriques-mutualises-pour-les-sciences-du-vivant/) led by the French Bioinformatics 
Institute ([IFB](https://www.france-bioinformatique.fr/)). It was initiated at the [IDRIS Open Hackathon](http://www.idris.fr/annonces/idris-gpu-hackathon-2023.html), part of the Open Hackathons program. 
The authors would like to acknowledge OpenACC-Standard.org for their support.