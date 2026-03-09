## Usage

Edit the `AFmassive_params.json` and/or `ColabFold_params.json` and/or `AlphaFold3_params.json` parameters file 
(see [file architecture](#tree)).  
Set first the [parameters of your run](https://github.com/GBLille/AFmassive?tab=readme-ov-file#running-afmassive) in the 
**AFM_run** section of the `AFmassive_params.json`, for instance:
```json
"AFM_run":
{
    "max_recycles": "20",
    "templates": "true",
    "dropout": "false",
    "dropout_structure_module": "false",
    "dropout_rates_filename": "",
    "stop_recycling_below": "0",
    "max_template_date": "2024-08-31",
    "min_score": "0",
    "max_score": "1",
    "fixed_seed": "",
    "db_preset": "full_dbs",
    "early_stop_tolerance": "0.5",
    "bfd_max_hits": "100000",
    "mgnify_max_hits": "501",
    "uniprot_max_hits": "50000",
    "uniref_max_hits": "10000"
}
```
or in the `ColabFold_params.json` file, for instance:
```json
"CF_run":
{
    "pair_strategy": "greedy",
    "use_dropout": "false",
    "num_recycle": "20",
    "recycle_early_stop_tolerance": "0.5",
    "stop_at_score": "100",
    "disable_cluster_profile": "false"
}
```
or in the `AlphaFold3_params.json` file, for instance
```json
"AF3_run":
{
    "fasta_chains": ["protein","protein"],
     "ligand": [
            {"ccdCodes": [""], "smiles": ""}
     ],
     "PTMs": [
            [{"type": "", "sequence": "", "positions": []}]
     ],
    "max_template_date": "2024-11-28",
    "num_diffusion_samples": "5",
    "unpairedMsa": "true",
    "pairedMsa": "true",
    "templates": "true"
}
```

**N.B.**: `"fasta_chains"` has to be filled before starting a run because no default value is provided. It specifies the 
type of each chain of the fasta file among `"protein"`, `"dna"` and `"rna"`. In this example, the two chains in the fasta 
file are proteins.

Then you can set the parameters of the **custom_params** section if necessary and the 
[plots section](#massivefold_plots-output-representation).

Activate the conda environment:
```bash
conda activate massivefold
```
>**N.B.**: on the Jean Zay cluster, simply load the `massivefold` >module. To be able to run on H100 or A100, uncomment the
>corresponding last lines of the `jobarray.slurm` header. Example >for H100:
>```bash
>module purge
>module load massivefold
>module load arch/h100
>```

Then launch MassiveFold:
```bash
massivefold run -s <SEQUENCE_PATH> -r <RUN_NAME> -p <NUMBER_OF_PREDICTIONS_PER_MODEL> -f <JSON_PARAMETERS_FILE> -t <TOOL>
```


Example for AFmassive:
```bash
massivefold run -s input/H1140.fasta -r afm_default -p 5 -f AFmassive_params.json
```
Example for ColabFold:
```bash
massivefold run -s input/H1140.fasta -r cf_default -p 5 -f ColabFold_params.json
```
Example for AlphaFold3:
```bash
massivefold run -s input/H1140.fasta -r af3_default -p 5 -f AlphaFold3_params.json
```

For more help and list of required and facultative parameters, run:
```bash
massivefold run -h
```
```text
usage: massivefold run [-h] -s SEQUENCE -r RUN_NAME -f PARAMETERS [-p PREDICTIONS_PER_MODEL] [-b BATCH_SIZE] [-j JOBID] [-o] [-c] [-C CALIBRATION_FROM] [-w WALL_TIME] [-m MSAS_PRECOMPUTED]
                       [-n TOP_N_MODEL] [-a] [--scheduler {auto,slurm,local}]

options:
  -h, --help            show this help message and exit

Required arguments:
  -s SEQUENCE, --sequence SEQUENCE
                        Path of the sequence(s) to infer, should be a 'fasta' file.
  -r RUN_NAME, --run RUN_NAME
                        Name chosen for the run to organize in outputs.
  -f PARAMETERS, --parameters PARAMETERS
                        Json file's path containing the parameters used for this run.

Optional arguments:
  -p PREDICTIONS_PER_MODEL, --predictions_per_model PREDICTIONS_PER_MODEL
                        Number of predictions (default: 5) computed for each neural network model. If used with -t AlphaFold3, -p is the number of seeds used. Each seed will have m samples
                        predicted. The number of sample set m is set in the AlphaFold3_params.json file. In total, with -p n, you will have m*n predictions computed.
  -b BATCH_SIZE, --batch_size BATCH_SIZE
                        Number of predictions per batch (default: 25). For AlphaFold3, it corresponds to number of seeds should not be higher than -p.
  -j JOBID, --jobid JOBID
                        Jobid of an alignment job to wait for inference, skips the alignments.
  -o, --only_msas       Only compute alignments, the first step of MassiveFold. Overwrite MSAs directory by forcing re-computation.
  -c, --calibrate       Calibrate --batch_size value. Searches from the previous runs for the same 'fasta' path given in --sequence and uses the longest prediction time found to compute
                        the maximal number of predictions per batch. This maximal number depends on the total time given by --wall_time.
  -C CALIBRATION_FROM, --calibration_from CALIBRATION_FROM
                        Path of a previous run to calibrate the batch size from (see --calibrate).
  -w WALL_TIME, --wall_time WALL_TIME
                        Total time in hour (default: 20) available for calibration computations.
  -m MSAS_PRECOMPUTED, --msas_precomputed MSAS_PRECOMPUTED
                        Path to directory that contains computed msas.
  -n TOP_N_MODEL, --top_n_model TOP_N_MODEL
                        Uses the n neural network models with best ranking confidence from this run's path.
  -a, --recompute_msas  Purges previous alignment step and recomputes msas.
  --scheduler {auto,slurm,local}
                        Scheduler selector (default: auto)
```

### Inference workflow

It launches MassiveFold with the same parameters introduced above but instead of running AFmassive, ColabFold or 
AlphaFold3 a single time, it divides it into multiple batches.

***N.B.***: with AFmassive and ColabFold, it is no more useful to specify if system to model is a monomer or a multimer 
as the detection is automatic now.

You can decide how the run will be divided by assigning `run_massivefold.sh` parameters *e.g.*:

```bash
massivefold run -s ./input/H1140.fasta -r 1005_preds -p 67 -b 25 -f AFmassive_params.json
```

The predictions are computed individually for each neural network (NN) model,  **-p** or **--predictions_per_model** 
allows to specify the number of predictions desired for each chosen model (number of seeds for AlphaFold3).  
These **--predictions_per_model** are then divided into batches with a fixed **-b** or **--batch_size** to optimize the 
run in parallel as each batch can be computed on a different GPU, if available.
The last batch of each NN model is generally smaller than the others to match the number of predictions fixed by 
**--predictions_per_model**.

***N.B.***: For each tool, the number of seeds is fixed by the **-p** parameter (for AFmassive and ColabFold, 
it means a different seed for each prediction; for AlphaFold3, it means a different seed for each samples set). 

For example, with **-b 25** and **-p 67** the predictions are divided into the following batches (separated runs), which 
are repeated for each NN model:

  1.  First batch: **--start_prediction=0** and **--end_prediction=24**
  2.  Second batch: **--start_prediction=25** and **--end_prediction=49**
  3.  Third batch: **--start_prediction=50** and **--end_prediction=67** 

By default (if **--models_to_use** is not assigned), all NN models are used: for protein complexes, 
15 models in total = 5 neural network models $\times$ 3 AlphaFold2 versions; for monomers, 5 
neural network models are used.

The prediction number per model can be adjusted, here with 67 per model and 15 models, it amounts to **1005 predictions 
in total divided into 45 batches**, these batches can therefore be run in parallel on a GPU cluster infrastructure.

The batch size can also be auto calibrated with the `-c` or `-C` parameters if at least one basic run has already been 
performed. The `-c` parameter will automatically search in the output folder that corresponds to the input sequence for 
the longest prediction duration. These options have to be coupled with the `-w` walltime parameter (it is advised to 
adapt this walltime value to the one of the job). For instance:

```bash
massivefold run -s ./input/H1140.fasta -r 1005_preds -p 67 -f AFmassive_params.json -c -w 10
```

***N.B.***: an interest to use `run_massivefold.sh` on a single server with a single GPU is to be able to run massive 
sampling for a structure in low priority, allowing other jobs with higher priority to be run in between.

### Parameters

#### Parameters in run_massivefold.sh

In addition to the parameters displayed with **-h** option, the json parameters file set with **-f** or **--parameters** 
should be organized like the `AFmassive_params.json` or `ColabFold_params.json` or `AlphaFold3_params.json` file.

#### Parameters in the json file

Each section of `AFmassive_params.json` or `ColabFold_params.json` or `AlphaFold3_params.json` is used for a different 
purpose.

The **massivefold** section designates the whole run parameters.  

```json
"massivefold": 
{
    "run_massivefold": "run_AFmassive.py",
    "run_massivefold_plots": "../massivefold/massivefold_plots.py",
    "data_dir": "$DSDIR/Alphafold-2024-04",
    "uniref_database": "",
    "jobfile_headers_dir": "./headers",
    "jobfile_templates_dir": "../massivefold/parallelization/templates",
    "output_dir": "./output",
    "logs_dir": "./log",
    "input_dir": "./input",
    "scripts_dir": "../massivefold/parallelization",
    "models_to_use": "",
    "pkl_format": "full"
}
```
The paths in the section are filled by `install.sh` but can be changed here if necessary. 
Headers (**jobfile_headers_dir**) are specified to setup the run, in order to give the parameters that are required to 
run the jobs on your cluster/server.
Build your own according to the [Jobfile's header building](#jobfiles-header-building) section.   
**models_to_use** is the list of NN models to use (for AFmassive and ColabFold). To select which NN models are used, 
separate them with a comma *e.g.*: "model_3_multimer_v1,model_3_multimer_v3", by default all are used  
**pkl_format**: how to manage pickle files (for AFmassive and ColabFold only)  
    - ‘full’ to keep the pickle files generated by the inference engine,  
    - ‘light’ to reduce its size by selecting main components, which are: number of recycles, PAE values, max PAE, 
plddt scores, ptm scores, iptm scores and ranking confidence values (stored in ./light_pkl directory)  
    - ‘none’ to remove them  
**uniref_database** is the parameter to fix the issue described in the [Troubleshooting](#troubleshooting) section (for AFmassive). 
You can change specifically uniref30 database path at this location. If not specified ("uniref_database": ""), the  
default uniref30 database path is used (same as AlphaFold2 configuration).  

- The **custom_params** section is relative to the personalized parameters that you want to add for your own cluster. 
For instance, for the Jean Zay GPU cluster:
```json
"custom_params": 
{
    "jeanzay_project": "<project>",
    "jeanzay_account": "<project>@v100",
    "jeanzay_gpu_with_memory": "v100-32g",
    "jeanzay_alignment_time": "10:00:00",
    "jeanzay_jobarray_time": "10:00:00"
}

```
As explained in [How to add a parameter](#how-to-add-a-parameter), these variables are substituted by their value when 
the jobfiles are created.

- For AFmassive, the **AFM_run** section gathers all the parameters used by MassiveFold for the run 
(see [AFmassive parameters](https://github.com/GBLille/AFmassive?tab=readme-ov-file#running-afmassive) 
section). All parameters except *--keep_pkl*, *--models_to_relax*, *--use_precomputed_msas*, *--alignment_only*, 
*--start_prediction*, *--end_prediction*, *--fasta_path* and *--output_dir* are exposed in this section.  
You can adapt the parameter values in function of your needs.  
The non exposed parameters mentioned before are set internally by the MassiveFold's pipeline or in the **massivefold** 
section (**models_to_use** and **pkl_format**).  

```json
"AFM_run":
{
    "max_recycles": "20",
    "templates": "true",
    "dropout": "false",
    "dropout_structure_module": "false",
    "dropout_rates_filename": "",
    "stop_recycling_below": "0",
    "max_template_date": "2024-08-31",
    "min_score": "0",
    "max_score": "1",
    "db_preset": "full_dbs",
    "early_stop_tolerance": "0.5",
    "bfd_max_hits": "100000",
    "mgnify_max_hits": "501",
    "uniprot_max_hits": "50000",
    "uniref_max_hits": "10000"
}
```
Lastly, the **plots** section is used for the MassiveFold plotting module. Only `DM_plddt_PAE` and `CF_PAEs` plots are 
available for AlphaFold3. The `recycles` plot is not available for AFmassive monomers. 
```json
"plots":
{
    "MF_plots_top_n_predictions": "10",
    "MF_plots_chosen_plots": "coverage,DM_plddt_PAE,CF_PAEs,score_distribution,recycles"
}
```

### Ranking

Three ranking files are generated for each run : `ranking_debug.json`, which contains the AlphaFold confidence score for each prediction, 
`ranking_ptm.json`, which contains the pTM score for each prediction, and, for complexes, `ranking_iptm.json`, which contains the ipTM score 
for each prediction. For AlphaFold3, a `confidences` subfolder contains detailed scores for each prediction, one json file per prediction.  

***N.B.***: For AFmassive and ColabFold, even if not directly comparable between AF2 versions because using weights coming from a different 
training, the scores are used to rank all the predictions for convenience, whatever the neural network version.

### Using ligands and modifications with AlphaFold3

In MassiveFold, ligands and post modifications are configured in the `AlphaFold3_params.json` file.  

For **ligands**, the `"ligand"` section has to be filled in with a CCD code **or** a SMILES code **or** a IUPAC code. In case of several, use several 
entries in the JSON as in the following example.  

For **modifications**, the `"modifications"` section has to be filled in. 
These are the available modifications as of yet:

| Name            | Chain type    | Target residue | Target base |
|-----------------|---------------|----------------|-------------|
| glycosylation   | protein       | N, S, T, K     | null        |
| phosphorylation | protein       | S, T           | null        |
| methylation     | protein & dna | R              | C           |
| hydroxylation   | protein       | P              | null        |
| acetylation     | protein       | K              | null        |
| cyclization     | protein       | E              | null        |

The `"modifications"` section contains as many entries (list) as the number of chains in the fasta file. The order of 
these chains is the same as in the fasta file and in the `"fasta_chains"` section.

For each modification, these two keys are required:
- `"type"` is the name of the modification (e.g.: `glycosylation`) 
- `"positions"` is a list of the positions on the (fasta) sequence where the modifications have to be linked

If the modification is a glycosylation, another key is needed:
- `"sequence"` is the sequence of the glycan in IUPAC code 

The following example shows 2 protein chains (the example H1140) with 4 ligands, a total of three glycosylations and 
one phosphorylation. The first protein chain is glycosylated once (residue 36) and phosphorylated once (residue 20), 
the second is glycosylated twice (same glycan on residues 74 and 84).

```json
"AF3_run":  {
        "fasta_chains": ["protein","protein"],
        "ligand": [
            {"ccdCodes": ["NAG"]},
            {"ccdCodes": ["KGM"]},
            {"smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"},
            {"IUPAC": "Gal(1-4)GlcNAc(1-2)Man(1-3)[Gal(1-4)GlcNAc(1-2)[Gal(1-4)GlcNAc(1-6)]Man(1-6)]Man(1-4)GlcNAc(1-4)[Fuc(1-6)]GlcNAc"}
        ],
        "modifications": [
            [
                    {
                        "type": "glycosylation",
                        "sequence": "Gal(1-4)GlcNAc(1-2)Man(1-3)[Gal(1-4)GlcNAc(1-2)[Gal(1-4)GlcNAc(1-6)]Man(1-6)]Man(1-4)GlcNAc(1-4)[Fuc(1-6)]GlcNAc",
                        "positions": [36]
                    },
                    {
                        "type": "phosphorylation",
                        "sequence": "",
                        "positions": [20]
                    }
            ],
            [
                    {
                        "type": "glycosylation",
                        "sequence": "Gal(1-4)GlcNAc(1-2)Man(1-3)[Gal(1-4)GlcNAc(1-2)[Gal(1-4)GlcNAc(1-6)]Man(1-6)]Man(1-4)GlcNAc(1-4)[Fuc(1-6)]GlcNAc",
                        "positions": [74,84]
                    }
            ]
        ],
        "max_template_date": "2024-11-28",
        "num_diffusion_samples": "5",
        "unpairedMsa": "true",
        "pairedMsa": "true",
        "templates": "true"
    }
```

### Relaxation

Because the relaxation takes time and resources to compute and that the MassiveFold process splits the predictions in 
many batches, the “use_gpu_relax” and “models_to_relax” parameters are set to “false” and “none” respectively. Indeed, 
if the relaxation is activated during the process, it will be run per batches, before the final ranking, resulting in 
relaxed structures that wouldn't necessarily be the best predictions. Instead, we recommend to use the `colabfold_relax` 
program provided in the `mf-colabfold` conda environment and developed by the ColabFold team, once all the predicted 
structures are produced and ranked. It allows to relax only selected PDB structures.  
For help, type:  
```bash
colabfold_relax -h
```

### Multiple runs gathering

We provide a `gather` command in the environment where MassiveFold is installed (also as a script at `src/massivefold/gather_runs.py`) that allows to collate the results of several runs. It 
gathers all the results and ranks them all. Run `gather -h` for help.

We also provide an `extract_scores.py` script that allows to extract the scores from pickle files and create rankings
(notably useful for interrupted runs). Run `python3 extract_scores.py -h` for help.

## Ligand screening
First, activate the conda environment:
```bash
conda activate massivefold
```
To launch a screening round, run:

```bash
massivefold screen -s <receptor_fasta_file> -l <ligand_list_csv> -f <AlphaFold3_params.json>
```

Run -h for help, as usual:
```bash
massivefold screen -h
```

### Format of the csv containing the ligands

This csv has 4 columns: `id`, `smiles`, `ccdCode` and `IUPAC`.  
Each row is a ligand to use for the screening round. "id" designates the name of the ligand, it can simply be the number of the ligand  in the list or a more detailed name for convenience.  
For the ligand sequence, use either `smiles`, `ccdCode` or `IUPAC` but not two at the same time, respectively in the SMILES format or the Chemical Compound Dictionnary code format (ccdCode) found at https://www.ebi.ac.uk/pdbe-srv/pdbechem/.

The format should be as follow:

| id       | smiles                   | ccdCode | IUPAC                  |
|----------|--------------------------|---------|------------------------|
| 1        |                          |         | Gal(1-4)GlcNAc(1-2)Man |
| 2        |                          | BGC     |                        |
| 3        |                          | OLA     |                        |
| aspirin  | CC(=O)OC1=CC=CC=C1C(=O)O |         |                        |
| glycerol |                          | GOL     |                        |

### Gathering the outputs

The predictions for each ligand will be located in a dedicated folder, itself located in an output folder 
named with the name of the csv file. All the results, notably the scores, can be gathered with the `gather` command (or script at `./src/massivefold/gather_runs.py`).  
See this [section](#multiple-runs-gathering).

## PPI screening 
First, activate the conda environment:
```bash
conda activate massivefold
```
To launch a PPI discovery round, run:

```bash
massivefold ppi --receptors <receptor_list> --ligands <ligand_list> -f <tool_params.json>
```
 
Run -h for help, as usual:
```bash
massivefold ppi -h
```
Additionally, a screening round (see [ligand screening](#ligand-screening)) can be applied on every PPI receptor-ligand pair by using the parameter `--context <ligand_csv>` (with the same format as [seen before](#format-of-the-csv-containing-the-ligands)).  
No need to change specify the `fasta_chains` field of the parameter file passed with `-f`, it will be automatically assigned in the PPI pipeline.

### Format of the csv for receptors and ligands

The `--receptors <receptor_list>` and `--ligands <ligand_list>` contain the absolute paths of the protein, DNA or RNA sequence(s) fasta files.  
This file is a CSV, with each line being a receptor (if passed to `--receptors`)  
or a ligand (if passed to `--ligands`) in this format:  
| protein              | DNA                 | RNA                 |
|----------------------|---------------------|---------------------|
| /path/to/prot1.fasta |                     |                     |
| /path/to/prot2.fasta |                     |                     |
| /path/to/prot3.fasta |                     |                     |
| /path/to/prot4.fasta |                     | /path/to/rna1.fasta |
| /path/to/prot5.fasta |                     |                     |
| /path/to/prot6.fasta |                     |                     |
|                      | /path/to/dna1.fasta |                     |
|                      |                     | /path/to/rna2.fasta |

