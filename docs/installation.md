## Installation

MassiveFold was developed to run massive sampling with [AFmassive](https://github.com/GBLille/AFmassive), 
[ColabFold](https://github.com/sokrypton/ColabFold) and [AlphaFold3](https://github.com/google-deepmind/alphafold3), and relies on them for its installation.

### Steps

<a id="retrieve-massivefold"></a>

1. **Retrieve MassiveFold**

```bash
# clone MassiveFold's repository
git clone https://github.com/GBLille/MassiveFold.git
```

For AFmassive runs, two additional installation steps are required to use MassiveFold:
- Download [sequence databases](https://github.com/GBLille/AFmassive?tab=readme-ov-file#sequence-databases)
- Retrieve the [neural network (NN) models parameters](https://github.com/GBLille/AFmassive?tab=readme-ov-file#alphafold-neural-network-model-parameters)

For ColabFold runs, two additional installation steps are required to use MassiveFold:
- Download [sequence databases](https://github.com/sokrypton/ColabFold?tab=readme-ov-file#generating-msas-for-large-scale-structurecomplex-predictions)
- Retrieve the [neural network (NN) models parameters](https://github.com/GBLille/AFmassive?tab=readme-ov-file#alphafold-neural-network-model-parameters)
and move them to a 'params' folder in the sequence databases folder (same parameters as AFmassive)

For AlphaFold3 runs, two additional installation steps are required to use MassiveFold:
- Download [sequence databases](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md#obtaining-genetic-databases)
- Retrieve the [neural network (NN) models parameters](https://github.com/google-deepmind/alphafold3?tab=readme-ov-file#obtaining-model-parameters)

<a id="install-massivefold"></a>

2. **Install MassiveFold**

We use an installation based on conda. The **install.sh** script we provide installs the conda environments using the 
`environment.yml`, `mf-afmassive.yml`, `mf-colabfold.yml` and `mf-alphafold3.yml` files. The first one is created for MassiveFold, 
the second one for AFmassive, the third for ColabFold and the last one for Alphafold3. It also creates the files 
architecture and set paths according to this architecture in the `AFmassive_params.json` and/or `ColabFold_params.json` 
and/or `AlphaFold3_params.json` parameters file.  

Help with:
```bash
cd MassiveFold
./install.sh -h
```

Installation with:
```bash
/install.sh [--only-envs] || --alphafold-db str --alphafold3-db str --colabfold-db str [--no-env]

./install -h for more details
  Options:
    --alphafold-db <str>: path to AlphaFold2 database
    --alphafold3-db <str>: path to AlphaFold3 database
    --colabfold-db <str>: path to ColabFold database
    --no-env: do not install the environments, only sets up the files and parameters.
      At least one of --alphafold-db or colabfold-db is required with this option.
    --only-envs: only install the environments (other arguments are not used)
```

This file tree displays the files' architecture after running `./install.sh`.

<a id="tree"></a> 
```txt
MassiveFold
├── install.sh
├── ...
├── examples
├── massivefold
└── massivefold_runs
    ├── AFmassive_params.json
    ├── ColabFold_params.json
    ├── AlphaFold3_params.json
    ├── headers/
        ├── example_header_alignment_jeanzay.slurm
        ├── example_header_jobarray_jeanzay.slurm
        └── example_header_post_treatment_jeanzay.slurm
    ├── input/
    ├── log/
    ├── output/
    └── run_massivefold.sh
```
The directory `massivefold_runs` is created, which contains:
- `AFmassive_params.json` to set the run parameters for AFmassive,
- `ColabFold_params.json` to set the run parameters for ColabFold,
- `AlphaFold3_params.json` to set the run parameters for AlphaFold3,
- `headers`' directory, containing the headers that must be created to use MassiveFold. Examples are given for the Jean 
Zay national CNRS French cluster (ready to use, see the [installation on Jean Zay](#install-on-jean-zay) to run 
MassiveFold directly on Jean Zay),
- `input` which contains the FASTA sequences,
- `log` with the logs of the MassiveFold runs (debug purposes), 
- `output` which contains the predictions, 
- `run_massivefold.sh` being the script to run [MassiveFold](#usage)

**On a GPU cluster:**  
- the **administrator** only needs to install the environments: 
```bash
cd MassiveFold
./install.sh --only-envs
```
- the **user** only needs to install the remaining files:
```bash
cd MassiveFold
./install.sh --no-env --alphafold-db <AF_DB_PATH> --colabfold-db <CF_DB_PATH> --alphafold3-db <AF3_DB_PATH>
```

***N.B.***: for <AF3_DB_PATH>, you need to use your own AlphaFold3 weights, while the database files can be shared. Therefore,
in a personal folder in your `home`, for instance in `~/af3_db` (that is your <AF3_DB_PATH>), you can copy your own weights, then 
in this folder, create symbolic links to the shared database files:

```bash
ls -s <SHARED_DB_PATH>/* ~/af3_db/
```

<a id="create-headers"></a>

3. **Create header files**  

Refer to [Jobfile's header](#jobfiles-header) for this installation step.

To run MassiveFold in parallel on your cluster/server, it is **required** to build custom jobfile headers for each step. 
They are three and should be named as follows: `{step}.slurm` (`alignment.slurm`, `jobarray.slurm` and 
`post_treatment.slurm`). The headers contain the parameters to give to SLURM for the jobs running (#SBATCH parameters). 
They have to be added in `MassiveFold/massivefold_runs/headers/` directory. Depending on your installation it can be 
another path, this path has to be set in the `AFmassive_params.json` and/or `ColabFold_params.json` and/or 
`AlphaFold3_params.json` as `jobfile_headers_dir` parameter.

Headers for Jean Zay cluster are provided as examples to follow (named `example_header_<step>_jeanzay.slurm`), to use 
them, rename each one following the previously mentioned naming convention.  

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

The jobfiles for each step are built by combining the jobfile header that you have to create in 
**MassiveFold/massivefold_runs/headers/** with the jobfile body in **massivefold/parallelization/templates/**.

Only the headers have to be adapted in function of your computation infrastructure. They contain the parameters to give to 
SLURM for the job running (#SBATCH parameters).
Each of the three headers (`alignment`, `jobarray` and `post treatment`) must be located in the **headers** directory 
(see [File architecture](#installation) section).

Their names should be identical to:
* **alignment.slurm**
* **jobarray.slurm**
* **post_treatment.slurm**

The templates work with the parameters provided in `AFmassive_params.json` and/or `ColabFold_params.json` and/or 
`AlphaFold3_params.json` files, given as a parameter to the **run_massivefold.sh** script.  
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
git clone https://github.com/GBLille/MassiveFold.git
./install.sh
```

To use AlphaFold3, copy your weights to the `~/af3_datadir` directory and run the following command:
```bash
ln -s $DSDIR/Alphafold3/* ~/af3_datadir/
```

The same [file architecture](#tree) is built, follow the [usage](#usage) section to use MassiveFold.

And specify the project you want to use in the `AFmassive_params.json` or `ColabFold_params.json` or `AlphaFold3_params.json`, 
replacing the `<project>` value by the 3-letters project name.  

***N.B.***: on Jean-Zay, AlphaFold3 runs only on A100 and H100. 

### Hardware recommendations

We recommend a 5 TB fast storage to host the sequence databases for AFmassive, ColabFold and AlphaFold3. The requirements in RAM 
depend on the length of the sequence(s) but 128GB should work in most cases. A GPU with at least 16 GB RAM is also recommended, 
knowing that more memory allows to model larger systems. 