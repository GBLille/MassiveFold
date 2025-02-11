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
    --tool <str>: (default: 'ColabFold') Use either AFmassive, ColabFold or ... in structure prediction for MassiveFold.
    --msas_precomputed <path>: path to directory that contains computed MSAs.
    --parameters <path>: json file's path containing the parameters used for this run.
    --predictions_per_model <int>: number of predictions computed for each neural network model.
    --batch_size <int>: (default: 25) number of predictions per batch, should not be higher than --predictions_per_model.
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

    // Treatment of the Multi-Sequence Alignment (MSA)
    def msa_results
    log.info("Running with: ${params.tool}")
    def seqFiles = files(params.sequence)


    // Check the Tool used and apply the appropriate workflow
    if (params.tool == 'ColabFold') {
        // Parse the JSON configuration
        def config = loadJson(params.config_tool)
        //println config.CF_run
        // FAIRE EN SORTE DE RECUPERE LES INFOS NECESSAIRES DANS UNE CHANNEL OU DANS PLUSIEURS ? 

        fastafile_unified = Unifier_colabfold(seqFiles)
        
        if (params.msas_precomputed) {
            // Use precomputed alignments
            log.info("Using precomputed MSAs: ${params.msas_precomputed}")
            msa_results = Channel
                .fromPath(params.msas_precomputed)
                .map { precomputed_msa -> tuple(precomputed_msa.baseName, precomputed_msa) }
        } else {
            // Run alignment process
            msa_results = Alignement_with_colabfold(fastafile_unified, params.database_dir, params.pair_strategy)
        }

        batch_json = Create_batchs(fastafile_unified, params.run, params.predictions_per_model, params.batch_size, params.config_tool, params.tool)

    // Process the batch JSON file to extract and structure batch information
    batchs = batch_json
        // Map each path in the batch_json channel to its corresponding data
        .map { path -> 
            // Convert the path to a File object
            File file = path.toFile()
            
            // Create a JSON parser to read the file
            def jsonSlurper = new JsonSlurper()
            
            // Parse the JSON file into a Groovy object
            def jsonData = jsonSlurper.parse(file)

            // Transform the JSON data into a list of tuples
            // Each entry in the JSON is converted to a tuple containing : key, start, end, model
            return jsonData.collect { key, value -> tuple(key, value.start, value.end, value.model) }
        }
        // Flatten the list of lists into a single list of tuples
        .flatten()
        // Group the flattened tuples into chunks of 4 elements each
        // This ensures each batch is represented as a single item in the channel
        .collate(4)


    batches_msa = batchs.combine(msa_results)

    res_prediction= Run_inference_colabfold(
            batches_msa,
            params.run,
            params.database_dir,
            params.num_recycle,
            params.recycle_early_stop_tolerance,
            params.use_dropout,
            params.stop_at_score,
            params.disable_cluster_profile
        )
    
    out = Standardize_output_colabfold(fastafile_unified, batch_json, res_prediction.collect())

    } 
    else if (params.tool == 'AFMassive') {
        //
        log.info("Do stuff")
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
            // To be continued 
        }
    } else {
        // Handle case where the tool is not recognized
        log.error("Unsupported tool: ${params.tool}")
        exit 1
    }


}

process Unifier_colabfold {
    input:
    path sequence

    output:
    path "input/converted_for_colabfold/${sequence.baseName}.fasta"

    script:
    """
    mkdir -p input
    cp $sequence input/.
    unifier.py --conversion input --to_convert input/${sequence.baseName}.fasta --tool ColabFold
    """
}

process Unifier_alphafold3 {
    input:
    path sequence

    output:
    path "${input_dir}/alphafold3_json_requests/${sequence.baseName}.fasta"

    script:
    """
    mkdir input
    cp $sequence input/.
    unifier.py --conversion input --to_convert input/${sequence.baseName}.fasta --tool AlphaFold3 --json_params AlphaFold3_params.json
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
    if [[ ${pair_strategy} == "greedy" ]]; then
        pairing_strategy=0
    elif [[ ${pair_strategy} == "complete" ]]; then
        pairing_strategy=1
    else
        echo "ValueError: --pair_strategy '${pair_strategy}'"
        exit 1
    fi

    colabfold_search $seqFile $data_dir ${seqFile.baseName}_msa --pairing_strategy \${pairing_strategy}
    """

    stub:
    """
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
        output_file="${seqFile.baseName}_\${sanitized_header}.a3m"
        # Write stub content
        echo "Stub content for \$line" > "\$output_file"
    done
    """
}


process RUN_alignment{
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
    if [[ ${pair_strategy} == "greedy" ]]; then
        pairing_strategy=0
    elif [[ ${pair_strategy} == "complete" ]]; then
        pairing_strategy=1
    else
        echo "ValueError: --pair_strategy '${pair_strategy}'"
        exit 1
    fi

    colabfold_search $seqFile $data_dir ${seqFile.baseName}_msa --pairing_strategy \${pairing_strategy}
    """

    stub:
    """
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
        output_file="${seqFile.baseName}_\${sanitized_header}.a3m"
        # Write stub content
        echo "Stub content for \$line" > "\$output_file"
    done
    """
}
process Create_batchs {
    input:
    path sequence
    val run
    val predictions_per_model
    val batch_size
    val config_tool
    val tool

    output:
    path "*.json"

    script:
    """
    batching.py --sequence_name=${sequence.baseName} \
                       --run_name=${run} \
                       --predictions_per_model=${predictions_per_model} \
                       --batch_size=${batch_size} \
                       --path_to_parameters=${config_tool} \
                       --tool ${tool}
    """
}


process Run_inference_colabfold {
    tag " $sequence_name | batch#$id_batch"
    label 'colabfold'

    input:
    tuple val(id_batch), val(batch_start), val(batch_end), val(batch_model), val(sequence_name) , path(msaFolder)
    val(run_name)
    val(data_dir)
    val(num_recycle)
    val(recycle_early_stop_tolerance)
    val(use_dropout)
    val(stop_at_score)
    val(disable_cluster_profile)

    output:
    //tuple val(id_batch), val(sequence_name), val(batch_start), val(batch_end), val(batch_model), path("res_*")
    path "res_*"

    script:
    """
    #export JAX_PLATFORMS=cpu
    random_seed=\$(shuf -i 0-1000000 -n 1)
    echo $batch_model
    num_models=\$(echo $batch_model | cut -d "_" -f 2)
    echo "using model \$num_models"
    num_version=\$(echo $batch_model | cut -d "v" -f 2)
    version_type="alphafold2_multimer_v\$num_version"
    model_type=\$version_type
    echo "from version version \$model_type"
    num_seeds=\$(($batch_end - $batch_start + 1))
    echo "predictions computed in the batch: \$num_seeds"
    echo "$msaFolder"
    colabfold_batch \
        ${msaFolder} \
        --data ${data_dir} \
        res_${sequence_name}_${run_name}_${id_batch} \
        --save-all \
        --random-seed \$random_seed \
        --num-seeds \$num_seeds \
        --model-type \$model_type \
        --model-order \$num_models \
        --num-models 1 \
        --num-recycle $num_recycle \
        --recycle-early-stop-tolerance $recycle_early_stop_tolerance \
        --stop-at-score $stop_at_score 
    """
        // $use_dropout \
        // $disable_cluster_profile
    
    // BOOL_use_dropout=$use_dropout
    // if $BOOl_disable_cluster_profile; then
    //     echo "Parameter --disable-cluster-profile set"
    //     disable_cluster_profile="--disable-cluster-profile"


    stub:
    """
    echo "Running in stub mode for batch $id_batch"
    # Simulate output files
    mkdir -p result_${sequence_name}_${run_name}_${id_batch}
    echo "Stub output for $sequence_name, batch $id_batch" > result_${sequence_name}_${run_name}_${id_batch}/prediction_stub.txt
    """
}


process Standardize_output_colabfold {
    publishDir "result/prediction/${seqFile.baseName}/"

    input: 
    path seqFile 
    path batchjson
    path colabfold_batch_files


    output: 
    path colabfold_batch_files

    script:
    """
    unifier.py --conversion output --to_convert . --batches_file $batchjson --tool ColabFold
    organize_outputs.py --batches_path .
    """
}

