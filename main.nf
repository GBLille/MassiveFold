nextflow.enable.dsl = 2
import groovy.json.JsonSlurper
import java.nio.file.Path
import java.nio.file.Paths

// Define the help message
def helpMessage = '''
Usage:
    nextflow main.nf --sequence <path> --run <name> --database_dir <path> [-options]
    
Required arguments:
    --sequence: path(s) of the sequence(s) to infer, should be a 'fasta' file or a list of files separated by commas.
    --run: name chosen for the run to organize outputs.
    --database_dir : path to the local directory where the tool database is located.

Optional arguments:
    --tool <str>: (default: 'ColabFold') Use either AFmassive, ColabFold or AlphaFold3 in structure prediction for MassiveFold.
    --msas_precomputed <path>: path to directory that contains computed MSAs.
    --parameters <path>: json file's path containing the parameters used for this run.
    --predictions_per_model <int>: (default: 1) number of predictions computed for each neural network model.
    --batch_size <int>: (default: 5) number of predictions per batch, should not be higher than --predictions_per_model.
    --pair_strategy <str>: (default: 'greedy') pairing strategy for MSA generation ('greedy' or 'complete').
    --num_recycle <int>: (default: 3) number of recycling iterations.
    --recycle_early_stop_tolerance <float>: (default: 0.5) early stop tolerance for recycling.
    --use_dropout <bool>: (default: false) use dropout during inference.
    --stop_at_score <float>: (default: 100) stop at this confidence score.
    --disable_cluster_profile <bool>: (default: false) disable cluster profile.
    --calibration_from <path>: path of a previous run to calibrate the batch size from (see --calibrate).
    --wall_time <int>: (default: 20) total time available for calibration computations, unit is hours.
    --top_n_models <int>: uses the n neural network models with best ranking confidence from this run's path.
    --only_msas: only compute alignments, the first step of MassiveFold.
    
Example in the environment of a virtual machine in the IFB-Biosphere cloud:
    nextflow main.nf --sequence examples/H1140.fasta --run test --database_dir ~/data/public/colabfold/ -profile docker -resume
'''

// Load the JSON file
def loadJson( fileName ) {
    def jsonSlurper = new JsonSlurper()
    return jsonSlurper.parse(new File(fileName))
}

workflow {
    // Display the help message if requested
    if (params.help || params.h) {
        log.info(helpMessage)
        exit 0
    }

    // Validate required parameters
    if (!params.sequence || !params.run || !params.database_dir) {
        log.error('Missing required parameters.')
        log.info(helpMessage)
        exit 1
    }

    // Debug information
    log.info("=== MassiveFold Pipeline ===")
    log.info("Tool: ${params.tool}")
    log.info("Sequence: ${params.sequence}")
    log.info("Run name: ${params.run}")
    log.info("Database directory: ${params.database_dir}")
    log.info("Pair strategy: ${params.pair_strategy}")
    log.info("Predictions per model: ${params.predictions_per_model}")
    log.info("Batch size: ${params.batch_size}")

    // Treatment of the Multi-Sequence Alignment (MSA)
    def msa_results
    def seqFiles = files(params.sequence)

    // Check the Tool used and apply the appropriate workflow
    if (params.tool == 'ColabFold') {
        // Parse the JSON configuration if it exists
        def config = [:]
        if (file(params.config_tool).exists()) {
            config = loadJson(params.config_tool)
            log.info("Loaded configuration from: ${params.config_tool}")
        } else {
            log.warn("Configuration file not found: ${params.config_tool}. Using default parameters.")
        }

        fastafile_unified = Unifier_colabfold(seqFiles)
        
        if (params.msas_precomputed) {
            // Use precomputed alignments
            log.info("Using precomputed MSAs: ${params.msas_precomputed}")
            msa_results = Channel
                .fromPath(params.msas_precomputed)
                .map { precomputed_msa -> tuple(precomputed_msa.baseName, precomputed_msa) }
        } else {
            // Run alignment process
            log.info("Running MSA alignment process")
            msa_results = Alignement_with_colabfold(fastafile_unified, params.database_dir, params.pair_strategy)
        }

        // Create batches
        batch_json = Create_batchs(
            fastafile_unified, 
            params.run, 
            params.predictions_per_model, 
            params.batch_size, 
            file(params.config_tool), 
            params.tool
        )

        // Parse batch JSON
        parsed_json = Parse_batch_json(batch_json)

        // Create batch channel
        batchs = parsed_json
            .splitCsv(sep: '\t')
            .map { row -> 
                def (id_batch, start, end, model) = row
                tuple(id_batch, start.toInteger(), end.toInteger(), model)
            }

        // Combine batches with MSA results
        batches_msa = batchs.combine(msa_results)

        // Run inference
        res_prediction = Run_inference_colabfold(
            batches_msa,
            params.run,
            params.num_recycle,
            params.recycle_early_stop_tolerance,
            params.use_dropout,
            params.stop_at_score,
            params.disable_cluster_profile
        )
        
        // Standardize output
        out = Standardize_output_colabfold(fastafile_unified, batch_json, res_prediction.collect())

    } 
    else if (params.tool == 'AFMassive') {
        log.error("AFMassive tool is not yet implemented")
        exit 1
    }
    else if (params.tool == 'AlphaFold3') {
        fastafile_unified = Unifier_AlphaFold3(seqFiles)
        if (params.msas_precomputed) {
            // Use precomputed alignments
            log.info("Using precomputed MSAs: ${params.msas_precomputed}")
            msa_results = Channel
                .fromPath(params.msas_precomputed)
                .map { precomputed_msa -> tuple(precomputed_msa.baseName, precomputed_msa) }
        } else {
            // Run alignment process
            msa_results = Alignement_with_colabfold(fastafile_unified, params.database_dir, params.pair_strategy)
            // TODO: Implement AlphaFold3 specific alignment
        }
        log.error("AlphaFold3 tool is not yet fully implemented")
        exit 1
    } else {
        // Handle case where the tool is not recognized
        log.error("Unsupported tool: ${params.tool}")
        log.error("Supported tools: ColabFold, AFMassive, AlphaFold3")
        exit 1
    }
}

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
    cp $sequence input/.
    unifier.py --conversion input --to_convert input/${sequence.baseName}.fasta --tool ColabFold
    """

    stub:
    """
    mkdir -p input/converted_for_colabfold
    cp $sequence input/converted_for_colabfold/${sequence.baseName}.fasta
    """
}

process Unifier_AlphaFold3 {
    label 'python_treatment'

    input:
    path sequence

    output:
    path "input/alphafold3_json_requests/${sequence.baseName}.fasta"

    script:
    """
    mkdir -p input/alphafold3_json_requests
    cp $sequence input/.
    unifier.py --conversion input --to_convert input/${sequence.baseName}.fasta --tool AlphaFold3 --json_params AlphaFold3_params.json
    """

    stub:
    """
    mkdir -p input/alphafold3_json_requests
    cp $sequence input/alphafold3_json_requests/${sequence.baseName}.fasta
    """
}

process Alignement_with_colabfold {
    tag "$seqFile.baseName"
    publishDir "result/${seqFile.baseName}/alignment"
    label 'colabfold'

    input:
    path(seqFile)
    path(data_dir)
    val(pair_strategy)

    output:
    tuple val(seqFile.baseName), path("${seqFile.baseName}_msa*")

    script:
    """

    echo "=== Starting ColabFold alignment ==="
    echo "Sequence file: $seqFile"
    echo "Data directory: $data_dir"
    echo "Pair strategy: ${pair_strategy}"
    ls -a > ls.txt

    # Check if data directory exists
    if [ ! -d "$data_dir" ]; then
        echo "ERROR: Database directory does not exist: $data_dir"
        exit 1
    fi
    
    # Set pairing strategy
    if [[ ${pair_strategy} == "greedy" ]]; then
        pairing_strategy=0
    elif [[ ${pair_strategy} == "complete" ]]; then
        pairing_strategy=1
    else
        echo "ValueError: --pair_strategy '${pair_strategy}' is not valid. Use 'greedy' or 'complete'"
        exit 1
    fi

    echo "Using pairing strategy: \$pairing_strategy"
    
    # Run ColabFold search
    colabfold_search $seqFile $data_dir ${seqFile.baseName}_msa --pairing_strategy \${pairing_strategy}
    
    echo "=== ColabFold alignment completed ==="
    """

    stub:
    """
    echo "=== Running in stub mode ==="
    # Create a directory for the MSA outputs
    msa_dir="${seqFile.baseName}_msa"
    mkdir -p "\$msa_dir"

    # Read the FASTA file and extract the headers
    grep '^>' $seqFile | while read -r line; do
        # Remove the leading '>'
        header=\$(echo "\$line" | sed 's/^>//')
        # Replace unsafe characters with underscores
        sanitized_header=\$(echo "\$header" | tr -c '[:alnum:]_' '_')
        # Create the output file name
        output_file="\$msa_dir/${seqFile.baseName}_\${sanitized_header}.a3m"
        # Write stub content
        echo "Stub MSA content for \$line" > "\$output_file"
        echo "Created stub file: \$output_file"
    done
    echo "=== Stub mode completed ==="
    """
}

process Create_batchs {
    label 'python_treatment'

    input:
    path sequence
    val run
    val predictions_per_model
    val batch_size
    path config_tool
    val tool

    output:
    path "*.json"

    script:
    """
    echo "=== Creating batches ==="
    echo "Sequence: ${sequence.baseName}"
    echo "Run: ${run}"
    echo "Predictions per model: ${predictions_per_model}"
    echo "Batch size: ${batch_size}"
    echo "Tool: ${tool}"
    
    batching.py --sequence_name=${sequence.baseName} \
                --run_name=${run} \
                --predictions_per_model=${predictions_per_model} \
                --batch_size=${batch_size} \
                --path_to_parameters=${config_tool} \
                --tool ${tool}
    
    echo "=== Batch creation completed ==="
    """

    stub:
    """
    echo "=== Creating stub batches ==="
    cat > ${sequence.baseName}_${run}_batches.json << 'EOF'
{
    "batch_1": {
        "start": 1,
        "end": ${batch_size},
        "model": "model_1_v3"
    }
}
EOF
    echo "=== Stub batch creation completed ==="
    """
}

process Run_inference_colabfold {
    tag "$sequence_name | batch#$id_batch"
    label 'colabfold'
    maxForks 2

    input:
    tuple val(id_batch), val(batch_start), val(batch_end), val(batch_model), val(sequence_name), path(msaFolder)
    val(run_name)
    
    val(num_recycle)
    val(recycle_early_stop_tolerance)
    val(use_dropout)
    val(stop_at_score)
    val(disable_cluster_profile)

    output:
    path "res_*"

    script:
    """
    echo "=== Starting ColabFold inference ==="
    echo "Batch ID: $id_batch"
    echo "Sequence: $sequence_name"
    echo "Batch range: $batch_start to $batch_end"
    echo "Model: $batch_model"
    echo "MSA folder: $msaFolder"
    
    # Generate random seed
    random_seed=\$(shuf -i 0-1000000 -n 1)
    echo "Random seed: \$random_seed"
    
    # Parse model information
    num_models=\$(echo $batch_model | cut -d "_" -f 2)
    echo "Model number: \$num_models"
    
    num_version=\$(echo $batch_model | cut -d "v" -f 2)
    version_type="alphafold2_multimer_v\$num_version"
    model_type=\$version_type
    echo "Model type: \$model_type"
    
    # Calculate number of seeds
    num_seeds=\$(($batch_end - $batch_start + 1))
    echo "Number of predictions: \$num_seeds"
    
    # Check MSA folder
    if [ ! -d "$msaFolder" ]; then
        echo "ERROR: MSA folder does not exist: $msaFolder"
        exit 1
    fi
    
    echo "MSA folder contents:"
    ls -la $msaFolder
    
    # Prepare additional flags
    additional_flags=""
    if [ "$use_dropout" = "true" ]; then
        additional_flags="\$additional_flags --use-dropout"
    fi
    if [ "$disable_cluster_profile" = "true" ]; then
        additional_flags="\$additional_flags --disable-cluster-profile"
    fi
    
    # Run ColabFold batch prediction
    colabfold_batch \
        ${msaFolder} \
        res_${sequence_name}_${run_name}_${id_batch} \
        --save-all \
        --random-seed \$random_seed \
        --num-seeds \$num_seeds \
        --model-type \$model_type \
        --model-order \$num_models \
        --num-models 1 \
        --num-recycle $num_recycle \
        --recycle-early-stop-tolerance $recycle_early_stop_tolerance \
        --stop-at-score $stop_at_score \
        \$additional_flags
    
    echo "=== ColabFold inference completed ==="
    """

    stub:
    """
    echo "=== Running inference in stub mode ==="
    echo "Batch ID: $id_batch"
    echo "Sequence: $sequence_name"
    
    # Create stub output directory
    output_dir="res_${sequence_name}_${run_name}_${id_batch}"
    mkdir -p "\$output_dir"
    
    # Create stub files
    echo "Stub prediction output for $sequence_name, batch $id_batch" > "\$output_dir/prediction_stub.txt"
    echo "Model: $batch_model" > "\$output_dir/model_info.txt"
    echo "Batch range: $batch_start to $batch_end" > "\$output_dir/batch_info.txt"
    
    echo "=== Stub inference completed ==="
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
    echo "=== Standardizing ColabFold outputs ==="
    echo "Sequence file: $seqFile"
    echo "Batch JSON: $batchjson"
    echo "Batch files: $colabfold_batch_files"
    
    # Create output directory
    mkdir -p organized_outputs
    
    # Run unifier
    unifier.py --conversion output --to_convert . --batches_file $batchjson --tool ColabFold
    
    # Organize outputs
    organize_outputs.py --batches_path .
    
    # Move organized outputs
    if [ -d "organized" ]; then
        mv organized/* organized_outputs/
    else
        # If organize_outputs.py doesn't create 'organized' directory, copy all results
        cp -r $colabfold_batch_files organized_outputs/ 2>/dev/null || true
    fi
    
    echo "=== Output standardization completed ==="
    """

    stub:
    """
    echo "=== Standardizing outputs in stub mode ==="
    mkdir -p organized_outputs
    echo "Stub standardized output for ${seqFile.baseName}" > organized_outputs/standardized_stub.txt
    echo "=== Stub standardization completed ==="
    """
}