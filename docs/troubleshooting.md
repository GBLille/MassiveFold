## Troubleshooting

### Uniref

While using AlphaFold2 or another AlphaFold2 based software like AFmassive, you can encounter a bug similar to this one in the 
msas generation:  

`WARNING: maximum number of residues 32763 exceeded in sequence UniRef100_A0A5E4B6R9_consensus`  
or  
`WARNING: maximum number of residues 32763 exceeded in sequence UniRef100_UPI000F443DA9 titin n=1 Tax=Lagenorhynchus 
obliquidens TaxID=90247 RepID=UPI000F443DA9`  
It makes the msas unusable and causes any inference on the sequence to crash. This was referenced in a 
[github issue](https://github.com/google-deepmind/alphafold/issues/810) and a [fix](https://github.com/google-deepmind/alphafold/issues/810#issuecomment-1666718050) was provided.  
According to this fix, download an updated version of the uniref30 database. To apply this modification to MassiveFold, 
set the `uniref_database` parameter in the AFmassive_params.json file to the updated database similarly to this:

```json
"massivefold": 
{
    "run_massivefold": "run_AFmassive.py",
    "run_massivefold_plots": "../massivefold/massivefold_plots.py",
    "data_dir": "/path/to/databases/AlphaFold/",
    "uniref_database": "/path/to/databases/AlphaFold/UniRef30_2023_02_hhsuite/UniRef30_2023_02",
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

### Usage without SLURM

MassiveFold can't run without SLURM. However, the `mf-afmassive-1.1.5`, `mf-colabfold-1.5.5` and `mf-alphafold-3.0.1` 
conda environments created at the installation allow to use respectively AFmassive, ColabFold and AlphaFold3 without parallelization. 
Their usage is detailed on their respective GitHub webpages.

### Alignment crashes with ColabFold

For a few sequences, the alignment step crashes for ColabFold. In this case, the web server can be called directly to get the alignments. 
For that, activate the ColabFold environment and run directly `colabfold_batch` with the `--msa-only` parameter. 
The alignment has to be located in a `msas_colabfold` folder of the output directory that corresponds to the fasta sequence to be able to run 
the structure inference with MassiveFold.  
In the `massivefold_runs` folder, run:

```commandline
conda activate mf-colabfold-1.5.5
colabfold_batch ./input/converted_for_colabfold/<sequence>.fasta ./output/<sequence>/msas_colabfold --msa-only
```

***N.B.***: The format of the fasta file to send to the ColabFold/MMseqs2 server has to be the MMseqs2 one, which means each sequence 
separated by ":".

### Using templates with ColabFold

MassiveFold is currently using the ColabFold 1.5.5 conda environment. It runs `colabfold_search` for the alignments and 
`colabfold_batch` for structure inference. Unfortunately, in this version, `colabfold_batch` doesn't allow to take a 
templates file as input, which means it would run the templates search for every inference job, which is not optimal.  
To be able to use templates with ColabFold, you can directly load the `mf-colabfold-1.5.5` environment and run 
`colabfold_batch` with the appropriate parameters, the `--save-all` parameter being mandatory (see `colabfold_batch -h`). 
With default parameters for alignments and templates research, it will query the MSA server.  
In case you would like to format the ColabFold output files to get the MassiveFold format, you need to put the ColabFold 
ouput files in a `batch_0` subfolder folder following the MassiveFold output file architecture :
```txt
output
└── <sequence name>
    └── <run name>
        └── batch_0 
```
Then, run the following scripts from the `massivefold_runs` folder, replacing `<sequence>` and `<run>` generic names 
with yours.
```commandline
python3 ../massivefold/parallelization/unifier.py \
    --to_convert output/<sequence>/<run>/batch_0/ \
    --conversion output_singular \
    --tool ColabFold

python3 ../massivefold/parallelization/organize_outputs.py \
    --batches_path output/<sequence>/<run>/
```

To generate plots, you can then use the `massivefold_plots.py` script as described [here](#massivefold_plots-output-representation).