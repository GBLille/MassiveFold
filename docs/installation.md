## Installation

MassiveFold was developed to run massive sampling with [AFmassive](https://github.com/GBLille/AFmassive), 
[ColabFold](https://github.com/sokrypton/ColabFold) and [AlphaFold3](https://github.com/google-deepmind/alphafold3), and relies on them for its installation.

### Steps

<a id="retrieve-massivefold"></a>

1. **Retrieve the databases and model parameters**

For AFmassive runs:
- Download [sequence databases](https://github.com/GBLille/AFmassive?tab=readme-ov-file#sequence-databases)
- Retrieve the [neural network (NN) models parameters](https://github.com/GBLille/AFmassive?tab=readme-ov-file#alphafold-neural-network-model-parameters)

For ColabFold runs:
- Download [sequence databases](https://github.com/sokrypton/ColabFold?tab=readme-ov-file#generating-msas-for-large-scale-structurecomplex-predictions)
- Retrieve the [neural network (NN) models parameters](https://github.com/GBLille/AFmassive?tab=readme-ov-file#alphafold-neural-network-model-parameters)
and move them to a 'params' folder in the sequence databases folder (same parameters as AFmassive)

For AlphaFold3 runs:
- Download [sequence databases](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md#obtaining-genetic-databases)
- Retrieve the [neural network (NN) models parameters](https://github.com/google-deepmind/alphafold3?tab=readme-ov-file#obtaining-model-parameters)

<a id="install-massivefold"></a>

2. **Install MassiveFold**

MassiveFold is a Python package, installed with pip in a conda environment created from the `environment.yml` file 
of the repository (Python and pip only):
```bash
conda env create -f environment.yml
conda activate massivefold
pip install massivefold
```
From a clone of the repository, `pip install .` installs that version instead. MassiveFold's environment only contains 
MassiveFold: the engines are installed in their own environments (see below). The jobs run MassiveFold with the Python 
of this environment, without activating it: it must be on a filesystem shared with the compute nodes.

Then `massivefold install` installs the structure prediction engines and creates the files architecture, with the 
paths of the databases set in the `AFmassive_params.json` and/or `ColabFold_params.json` and/or `AlphaFold3_params.json` 
parameters files:
```bash
massivefold install --alphafold-db <AF_DB_PATH> --colabfold-db <CF_DB_PATH> --alphafold3-db <AF3_DB_PATH> [--install-path massivefold_runs]
```
Only the engines of the databases given are installed.

```txt
massivefold install -h
  --alphafold-db <str>: path to AlphaFold2 database
  --alphafold3-db <str>: path to AlphaFold3 database
  --colabfold-db <str>: path to ColabFold database
  --install-path <str>: where to create the MassiveFold file architecture, engines included (default: massivefold_runs)
  --no-env: do not install the environments, only files and parameters: the engines are run from the PATH.
  --only-envs: only install the environments of the engines (of the given databases, or all of them)
  --jeanzay: install for the Jean Zay cluster (detected automatically on its login nodes and in its jobs)
  --ignore-conda-envs: install all the engines' environments, even those already found in conda
```

**The engines** (AFmassive, ColabFold and AlphaFold3) each run in their own environment, as they need different and 
incompatible library versions. They are installed with [pixi](https://pixi.sh), from the recipes 
(`massivefold/engines/pixi.toml`, locked by `pixi.lock`). The jobs run each engine through a launcher, 
`<engines-dir>/bin/run_in_<engine>`, recorded as `engine_prefix` in the parameters files. To use engines installed 
another way, set `engine_prefix` to the command that runs them (for instance a container `apptainer exec --nv <image.sif>`), 
or leave it empty if they are in the `PATH` of the jobs (`--no-env`).

An engine already installed in conda, in an environment named (`mf-afmassive-1.1.11`, `mf-colabfold-1.6.3`, 
`mf-alphafold3-1.2.0` for this MassiveFold version), is used instead of being installed. Use `--ignore-conda-envs` 
to install all the engines with pixi anyway.

The jobs run each engine through a launcher, `<install-path>/engines/bin/run_in_<engine>`, recorded as `engine_prefix` in the 
parameters files. It runs the command in the engine's environment, without pixi or network access on the compute 
nodes. To use engines installed another way, set `engine_prefix` to the command that runs them (for instance a container 
`apptainer exec --nv <image.sif>`), or leave it empty if they are in the `PATH` of the jobs (`--no-env`).

This file tree displays the files' architecture after running `massivefold install`.

<a id="tree"></a> 
```txt
massivefold_runs
├── AFmassive_params.json
├── ColabFold_params.json
├── AlphaFold3_params.json
├── multirun_setup.json
├── engines/
|    ├── bin/
|    |    ├── run_in_afmassive
|    |    ├── run_in_colabfold
|    |    └── run_in_alphafold3
|    ├── pixi.toml
|    └── pixi.lock
├── headers/
|    ├── alignment.slurm
|    ├── jobarray.slurm
|    ├── post_treatment.slurm
|    ├── example_header_alignment_jeanzay.slurm
|    ├── example_header_jobarray_jeanzay.slurm
|    └── example_header_post_treatment_jeanzay.slurm
├── input/
|    └── H1140.fasta
├── log/
└── output/
```
The directory `massivefold_runs` is created, which contains:
- `AFmassive_params.json` to set the run parameters for AFmassive,
- `ColabFold_params.json` to set the run parameters for ColabFold,
- `AlphaFold3_params.json` to set the run parameters for AlphaFold3,
- `engines` with the environments of the engines,
- `headers`' directory, containing the jobfile headers (see [Create header files](#create-headers)),
- `input` which contains the FASTA sequences,
- `log` with the logs of the MassiveFold runs (debug purposes), 
- `output` which contains the predictions

The database paths and the engines' launchers are also kept in `~/.config/massivefold/site_defaults/`, so that a later 
`massivefold install` of another workspace reuses them.

**On a GPU cluster:**  
- the **administrator** only needs to install the engines once, in a shared directory: 
```bash
massivefold install --only-envs --install-path <SHARED_DIR>
```
- the **user** installs MassiveFold and only creates the remaining files, then sets `engine_prefix` in the parameters 
files to the shared launchers (`<SHARED_DIR>/engines/bin/run_in_<engine>`):
```bash
massivefold install --no-env --alphafold-db <AF_DB_PATH> --colabfold-db <CF_DB_PATH> --alphafold3-db <AF3_DB_PATH>
```

***N.B.***: for <AF3_DB_PATH>, you need to use your own AlphaFold3 weights, while the database files can be shared. Therefore,
in a personal folder in your `home`, for instance in `~/af3_db` (that is your <AF3_DB_PATH>), you can copy your own weights, then 
in this folder, create symbolic links to the shared database files:

```bash
ln -s <SHARED_DB_PATH>/* ~/af3_db/
```

<a id="create-headers"></a>

3. **Adapt the header files**  

Refer to [Jobfile's header](#jobfiles-header) for this installation step.

To run MassiveFold in parallel on your cluster/server, jobfile headers are needed for each step. They are three and 
named as follows: `{step}.slurm` (`alignment.slurm`, `jobarray.slurm` and `post_treatment.slurm`). The headers contain 
the parameters to give to SLURM for the jobs running (#SBATCH parameters). They are in the 
`massivefold_runs/headers/` directory. Depending on your installation it can be another path, this path has to be set in 
the `AFmassive_params.json` and/or `ColabFold_params.json` and/or `AlphaFold3_params.json` as `jobfile_headers_dir` 
parameter.

`massivefold install` creates generic headers as a basis to adapt for your cluster (nodes, RAM, partitions, time limits...). 
Headers for Jean Zay cluster are provided as examples to follow (named `example_header_<step>_jeanzay.slurm`).
<a id="set-custom-params"></a>

4. **Set custom parameters**

Each cluster has its own specifications in parameterizing job files. For flexibility needs, you can add your custom 
parameters in your headers, and then in the `AFmassive_params.json` file and/or `ColabFold_params.json` file and/or 
`AlphaFold3_params.json` so that you can dynamically change their values in the json file.  

To illustrate these "special needs", here is an example of parameters that can be used on the French national Jean Zay 
cluster to specify GPU type, time limits or the project on which the hours are used:

Go to `AFmassive_params.json`/`ColabFold_params.json`/`AlphaFold3_params.json` location and modify it:
```bash
cd MassiveFold/massivefold_runs

"custom_params":
{
    "jeanzay_gpu": "v100",
    "jeanzay_project": "<project>",
    "jeanzay_account": "<project>@v100",
    "jeanzay_gpu_with_memory": "v100-32g",
    "jeanzay_alignment_time": "05:00:00",
    "jeanzay_jobarray_time": "15:00:00"
}
```
And specify them in the jobfile headers (such as here for `MassiveFold/headers/jobarray.slurm`) 
```
#SBATCH --account=$jeanzay_account
#SBATCH -C $jeanzay_gpu_with_memory
#SBATCH --time=$jeanzay_jobarray_time
```
### Jobfile's header

#### Building

The jobfiles for each step are built by combining the jobfile header in **massivefold_runs/headers/** with the 
jobfile body in **massivefold/parallelization/templates/** (in the installed package).

Only the headers have to be adapted in function of your computation infrastructure. They contain the parameters to give to 
SLURM for the job running (#SBATCH parameters).
Each of the three headers (`alignment`, `jobarray` and `post treatment`) must be located in the **headers** directory 
(see [File architecture](#installation) section).

Their names should be identical to:
* **alignment.slurm**
* **jobarray.slurm**
* **post_treatment.slurm**

The templates work with the parameters provided in `AFmassive_params.json` and/or `ColabFold_params.json` and/or 
`AlphaFold3_params.json` files, given as a parameter to the `massivefold` command.  
These parameters are substituted in the template job files thanks to the python library [string.Template](https://docs.python.org/3.8/library/string.html#template-strings).  
Refer to [How to add a parameter](#how-to-add-a-parameter) for parameters substitution.

- **Requirement:** In the jobarray's jobfile header (*massivefold_runs/headers/jobarray.slurm*) should be stated that 
it is a job array and the number of tasks in it has to be given. The task number argument is substituted with the 
*$substitute_batch_number* parameter.  
It should be expressed as:
```
#SBATCH --array=0-$substitute_batch_number
```
For example, if there are 45 batches, with 1 batch per task of the job array, the substituted expression will be:
```
#SBATCH --array=0-44
```
- Add these lines too in the headers, it is necessary to store MassiveFold's log:

In **alignment.slurm**:
```
#SBATCH --error=${logs_dir}/${sequence_name}/${run_name}/alignment.log
#SBATCH --output=${logs_dir}/${sequence_name}/${run_name}/alignment.log
```
In **jobarray.slurm**:

```
#SBATCH --error=${logs_dir}/${sequence_name}/${run_name}/jobarray_%a.log
#SBATCH --output=${logs_dir}/${sequence_name}/${run_name}/jobarray_%a.log
```
In **post_treatment.slurm**:
```
#SBATCH --output=${logs_dir}/${sequence_name}/${run_name}/post_treatment.log
#SBATCH --error=${logs_dir}/${sequence_name}/${run_name}/post_treatment.log
```
We provide headers for the Jean Zay French CNRS national GPU cluster ([IDRIS](http://www.idris.fr/),) 
that can also be used as examples for your own infrastructure.

#### How to add a parameter

- Add **\$new_parameter** or **\$\{new_parameter\}** in the template's header where you want its value to be set and 
in the "custom_params" section of `AFmassive_params.json` and/or `ColabFold_params.json` and/or `AlphaFold3_params.json` 
where its value can be specified and modified conveniently for each run.

**Example** in the json parameters file for Jean Zay headers:
```json
"custom_params":
{
    "jeanzay_account": "project@v100",
    "jeanzay_gpu_with_memory": "v100-32g",
    "jeanzay_jobarray_time": "10:00:00"
}
```
Where "project" is your 3 letter project with allocated hours on Jean Zay.

- These parameters will be substituted in the header where the parameter keys are located:

```
#SBATCH --account=$jeanzay_account

#SBATCH --error=${logs_dir}/${sequence_name}/${run_name}/jobarray_%a.log
#SBATCH --output=${logs_dir}/${sequence_name}/${run_name}/jobarray_%a.log

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --gpus-per-node=1
#SBATCH --array=0-$substitute_batch_number
#SBATCH --time=$jeanzay_jobarray_time
##SBATCH --qos=qos_gpu-t4               # Uncomment for job requiring more than 20h (max 16 GPUs)
#SBATCH -C $jeanzay_gpu_with_memory     # GPU type+memory
```
- Never use single \$ symbol for other uses than parameter/value substitution from the json file.\
To use $ inside the template files (bash variables or other uses), use instead $$ as an escape following 
[string.Template](https://docs.python.org/3.8/library/string.html#template-strings) documentation.

### Installation on Jean Zay

To use it on Jean Zay, the only installation steps are:
```bash
module load miniforge/26.1.1-3
conda env create -f environment.yml -p $WORK/conda/massivefold
conda activate $WORK/conda/massivefold
pip install massivefold
massivefold install [--install-path str]
```
On Jean Zay, the prebuilt headers and parameters files are used, all the engines are installed and the AlphaFold3 
databases of `$DSDIR/Alphafold3` are linked in `~/af3_datadir`. 
To use AlphaFold3, **copy your weights to the `~/af3_datadir` directory**.

The same [file architecture](#tree) is built, follow the [usage](#usage) section to use MassiveFold.
As the engines' environments are large, use an `--install-path` in your `$WORK` rather than in your `$HOME`.
And specify the project you want to use in the `AFmassive_params.json` or `ColabFold_params.json` or `AlphaFold3_params.json`, 
replacing the `<project>` value by the 3-letters project name.  

***N.B.***: on Jean-Zay, AlphaFold3 runs only on A100 and H100. 

### Hardware recommendations

We recommend a 5 TB fast storage to host the sequence databases for AFmassive, ColabFold and AlphaFold3. The requirements in RAM 
depend on the length of the sequence(s) but 128GB should work in most cases. A GPU with at least 16 GB RAM is also recommended, 
knowing that more memory allows to model larger systems. 
