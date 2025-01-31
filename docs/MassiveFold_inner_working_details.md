# Defining variables:
  - {tool}_params.json: json file containing all params given as an input (-f). The name can be anything.
  - {scripts_dir}: defined in the {tool}_params.json in "massivefold" section.
  - {input}: defined in the {tool}_params.json in "massivefold" section.
  - {output}: defined in the {tool}_params.json in "massivefold" section.
  - {log_dir}: defined in the {tool}_params.json in "massivefold" section. Move every files containing logs inside.
  - {data_dir}: defined in the {tool}_params.json in "massivefold" section. Contain databases and NN weights.
  - {models_to_use}: defined in the {tool}_params.json in "massivefold" section. Names of the AF2 NN models to use for the predictions.
  - {sequence}: name of the fasta sequence given as an input (-s), extracted as its basename without the ".fasta".
  - {run}: name of the run chosen by the user, given as an input (-r)..
  - {tool}: name of the tool to use given as an input (-t), Limited to values "AFmassive", "ColabFold", "AlphaFold3" for now.
  - {predictions_per_model}: given as an input (-p), number of prediction per NN model (for AF2 based model) or diffusion seed (AF3).
  - {batch_size}: given as an input (-b), number of predictions (-p) to fit in a single batch.
  - {concat_sequences}: arbitrary name given by ColabFold for the a3m msas file. Irrelevant in the process as the filename is used nowhere.
  - {batch_[i]_start}, {batch_[i]_end}, {batch_[i]_model}, value of the batch i entries in the {sequence}\_{run}_batches.json file., 

# Describing each step from the inner working diagram
- unifier.py to use standard fasta format with tools that demand specific format  
  A. ```bash python ./{scripts_dir}/unifier.py  --conversion input --to_convert ./{input}/{sequence}.fasta --tool ColabFold```  
    Input: {input}/{sequence}.fasta  
    Output: {input}/converted_for_colabfold/{sequence}.fasta  

  B. ```bash python ./{scripts_dir}/unifier.py  --conversion input --to_convert ./{input}/{sequence}.fasta --tool AlphaFold3 --json_params AlphaFold3_params.json```  
    Input: {sequence}.fasta and {AlphaFold3_params}.json  
    Output: {input}/alphafold3_json_requests/{sequence}.fasta  

- Alignment on CPU:  
  C. ```bash run_AFmassive.py --alignments_only=true --fasta_paths={input}/{sequence}.fasta --output_dir={output}/```  
    Input: {input}/{sequence}.fasta  
    Output: {output}/{sequence}/msas/chain_id_map.json  
    N.B. : Output is not limited to these file, there is also one directory per chain in {output}/{sequence}/.

  D. ```bash colabfold_search {input}/converted_for_colabfold/{sequence}.fasta {data_dir} {output}/{sequence}/msas_colabfold/```  
    Input: {input}/converted_for_colabfold/{sequence}.fasta  
    Output: {output}/{sequence}/msas_colabfold/{concat_sequences}.a3m  

  E. ```bash run_alphafold.py --norun_inference --json_path {input}/alphafold3_json_requests/{sequence}.json --output_dir {output}/{sequence}```  
    Input: {input}/alphafold3_json_requests/{sequence}.fasta  
    Output: {output}/{sequence}/msas_alphafold3/msas_alphafold3_data.json  

- Batch creation (all tools) and batch setup (AlphaFold3)  
  F. ```bash python {scripts_dir}/batching.py  
          --sequence_name={sequence}  
          --run_name={run}  
          --predictions_per_model={predictions_per_model}  
          --batch_size={batch_size}  
          --models_to_use={models_to_use}  
          --path_to_parameters={tool}_params.json  
          --tool {tool}```  
    Input: {tool}_params.json, ...(others are described in the command above.)  
    Output: {sequence}\_{run}_batches.json  
    N.B. : After being used for the following jobs as an input, this json file containing batches information is moved to the {log_dir}/{sequence}/{run}/ directory.  

  G. ```bash python {scripts_dir}/unifier.py  
          --conversion input_inference  
          --to_convert {output}/{sequence}/msas_alphafold3/msas_alphafold3_data.json  
          --json_params {tool}_params.json  
          --batches_file {sequence}_{run}_batches.json  
          --tool AlphaFold3```  
    Input: Described in the command above.  
    Output: {output}/{sequence}/{run}/af3_batch_[i].json (With i going from 0 to n, n being the number of batches minus one)  

- Inference on GPU  
  H. ```bash run_AFmassive.py --use_precomputed_msa=trues --start_prediction={batch_[i]_start} --end_prediction={batch_[i]_end} --models_to_use={batch_[i]_model}```  
    Input: Batches information that are stored in the {sequence}\_{run}\_batches.json. They can be easily retrieved with the {script_dir}/get_batch.py script.  
    Output: {output}/{sequence}/{run}/batch_[i]/ (With i going from 0 to n, n being the number of batches minus one)  

  I. ```bash colabfold_batch  
          {output}/{sequence}/msas_colabfold/  
          {output}/batch_[i]  
          --save-all  
          --num-models 1  
          --num-seeds ({batch_[i]_end} - {batch_[i]_start} + 1)```  
    Input: Batches information that are stored in the {sequence}\_{run}\_batches.json. They can be easily retrieved with the {script_dir}/get_batch.py script.  
    Output: {output}/{sequence}/{run}/batch_[i]/ (With i going from 0 to n, n being the number of batches minus one)  
    N.B: The output in each 'batch_[i]' subdirectory is formatted as any ColabFold output  

  J. ```bash run_alphafold.py --norun_data_pipeline --json_path {output}/{sequence}/{run}/af3_batch_[i].json```  
    Input: {output}/{sequence}/{run}/af3_batch_[i].json (With i going from 0 to n, n being the number of batches minus one)  
    Output: {output}/{sequence}/{run}/batch_[i]/ (With i going from 0 to n, n being the number of batches minus one)  
    N.B: The output in each 'batch_[i]' subdirectory is formatted as any AlphaFold3 output  

- Standardize the outputs from original tool format to AF2 standard format  
  K ```bash python {scripts_dir}/unifier.py --conversion output--to_convert {output}/{sequence}/{run}/ --batches_file {sequence}_{run}_batches.json --tool ColabFold```  
    Input: Batch directory, each batch containing ColabFold output  
    Output: ColabFold output in each batch subdirectories formatted as the standard AF2 output  

  L. ```bash python {scripts_dir}/unifier.py --conversion output--to_convert {output}/{sequence}/{run}/ --batches_file {sequence}_{run}_batches.json --tool AlphaFold3```  
    Input: Batch directory, each batch containing AlphaFold3 output  
    Output: AlphaFold3 output in each batch subdirectories formatted as the standard AF2 output  

- Organize the batches outputs to makes them united in one output  
  M. ```bash python {scripts_dir}/organize_outputs.py --batches_path {output}/{sequence}/{run}/```  

- Lighten the output size  
  N. ```bash python {scripts_dir}/lighten_output.py {output}/{sequence}/{run}/```  
