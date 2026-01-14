nextflow.enable.dsl = 2

/*
 * =====================================
 *   Processes à l’outil ColabFold
 * =====================================
 */

process Parse_batch_json {
    label 'python_treatment'

    input:
    path batch_json_file

    output:
    path "parsed_batches.tsv"

    script:
    """
    python3 -c "
import json

with open('${batch_json_file}', 'r') as f:
    data = json.load(f)

with open('parsed_batches.tsv', 'w') as out:
    for key, val in data.items():
        start = val.get('start')
        end = val.get('end')
        model = val.get('model')
        if start is not None and end is not None and model is not None:
            out.write(f'{key}\\t{start}\\t{end}\\t{model}\\n')
        else:
            print(f'Skipping incomplete entry for key: {key}')
    "
    """
}

process Unifier_colabfold {
    label 'python_treatment'

    input:
    path sequence

    output:
    path "input/converted_for_colabfold/${sequence.baseName}.fasta"

    script:
    """
    mkdir -p input/converted_for_colabfold
    cp ${sequence} input/.
    unifier.py --conversion input --to_convert input/${sequence.baseName}.fasta --tool ColabFold
    """

    stub:
    """
    mkdir -p input/converted_for_colabfold
    cp ${sequence} input/converted_for_colabfold/${sequence.baseName}.fasta
    """
}

process Alignement_with_colabfold {
    tag "$seqFile.baseName"
    publishDir "result/${seqFile.baseName}/alignment"
    label 'colabfold'

    input:
    path(seqFile)
    val(data_dir)
    val(pair_strategy)

    output:
    tuple val(seqFile.baseName), path("${seqFile.baseName}_msa")

    script:
    """
    echo "=== Starting ColabFold alignment ==="
    echo "Sequence file: ${seqFile}"
    echo "Data directory: ${data_dir}"
    echo "Pair strategy: ${pair_strategy}"

    if [ ! -d "${data_dir}" ]; then
        echo "ERROR: Database directory does not exist: ${data_dir}"
        exit 1
    fi
    
    if [[ "${pair_strategy}" == "greedy" ]]; then
        pairing_strategy=0
    elif [[ "${pair_strategy}" == "complete" ]]; then
        pairing_strategy=1
    else
        echo "ValueError: --pair_strategy '${pair_strategy}' is not valid. Use 'greedy' or 'complete'"
        exit 1
    fi

    colabfold_search ${seqFile} ${data_dir} ${seqFile.baseName}_msa --pairing_strategy \${pairing_strategy}
    echo "=== ColabFold alignment completed ==="
    """

    stub:
    """
    echo "=== Running alignment in stub mode ==="
    msa_dir="${seqFile.baseName}_msa"
    mkdir -p "\${msa_dir}"

    grep '^>' ${seqFile} | while read -r line; do
        header=\$(echo "\${line}" | sed 's/^>//')
        sanitized_header=\$(echo "\${header}" | tr -c '[:alnum:]_' '_')
        output_file="\${msa_dir}/${seqFile.baseName}_\${sanitized_header}.a3m"
        echo "Stub MSA content for \${line}" > "\${output_file}"
    done
    """
}

process Create_batchs {
    label 'python_treatment'

    input:
    path sequence
    val run
    val predictions_per_model
    val batch_size
    val tool

    output:
    path "*.json"

    script:
    def toolCode = tool.toLowerCase()
    def configMap = [
        massivefold: params.massivefold,
        CF_run: params.CF_run,
        AFM_run: params.AFM_run,
        plots: params.plots
    ]
    def configJson = groovy.json.JsonOutput.toJson(configMap)
    
    """
    cat > temp_config.json << 'EOF'
${configJson}
EOF

    batching.py --sequence_name=${sequence.baseName} \
                --run_name=${run} \
                --predictions_per_model=${predictions_per_model} \
                --batch_size=${batch_size} \
                --path_to_parameters=temp_config.json \
                --tool ${tool}
    """
}

process Run_inference_colabfold {
    tag "$sequence_name | batch#$id_batch"
    label 'colabfold'
    maxForks 2

    input:
    tuple val(id_batch), val(batch_start), val(batch_end), val(batch_model), val(sequence_name), path(msaFolder)
    val run_name
    val num_recycle
    val recycle_early_stop_tolerance
    val use_dropout
    val stop_at_score
    val disable_cluster_profile

    output:
    path "res_*"

    script:
    """
    if command -v nvidia-smi &> /dev/null && nvidia-smi -L | grep -q "GPU"; then
        export CUDA_VISIBLE_DEVICES=0
    else
        export CUDA_VISIBLE_DEVICES=""
        export JAX_PLATFORMS=cpu
    fi

    random_seed=\$(shuf -i 0-1000000 -n 1)
    num_models=\$(echo ${batch_model} | cut -d "_" -f 2)
    num_version=\$(echo ${batch_model} | cut -d "v" -f 2)
    model_type="alphafold2_multimer_v\${num_version}"
    num_seeds=\$((${batch_end} - ${batch_start} + 1))

    additional_flags=""
    if [ "${use_dropout}" = "true" ]; then
        additional_flags="\${additional_flags} --use-dropout"
    fi
    if [ "${disable_cluster_profile}" = "true" ]; then
        additional_flags="\${additional_flags} --disable-cluster-profile"
    fi

    colabfold_batch \
        ${msaFolder} \
        res_${sequence_name}_${run_name}_${id_batch} \
        --save-all \
        --random-seed \${random_seed} \
        --num-seeds \${num_seeds} \
        --model-type \${model_type} \
        --model-order \${num_models} \
        --num-models 1 \
        --num-recycle ${num_recycle} \
        --recycle-early-stop-tolerance ${recycle_early_stop_tolerance} \
        --stop-at-score ${stop_at_score} \
        \${additional_flags}
    """

    stub:
    """
    mkdir -p "res_${sequence_name}_${run_name}_${id_batch}"
    echo "Stub prediction output" > "res_${sequence_name}_${run_name}_${id_batch}/prediction_stub.txt"
    """
}

process Standardize_output_colabfold {
    publishDir "result/prediction/${seqFile.baseName}/"
    label 'python_treatment'

    input: 
    path seqFile 
    path batchjson
    path colabfold_batch_files

    output: 
    path "organized_outputs/*"

    script:
    """
    mkdir -p organized_outputs
    unifier.py --conversion output --to_convert . --batches_file ${batchjson} --tool ColabFold
    organize_outputs.py --batches_path .

    if [ -d "organized" ]; then
        mv organized/* organized_outputs/
    else
        cp -r res_* organized_outputs/ 2>/dev/null || true
    fi
    """

}
