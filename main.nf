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
    --config_json <path>: (default: './config.json') path to JSON configuration file.
    --msas_precomputed <path>: path to directory that contains computed MSAs.
    --predictions_per_model <int>: (default: 1) number of predictions computed for each neural network model.
    --batch_size <int>: (default: 5) number of predictions per batch, should not be higher than --predictions_per_model.
    --pair_strategy <str>: (default: 'greedy') pairing strategy for MSA generation ('greedy' or 'complete').
    --num_recycle <int>: (default: 3) number of recycling iterations.
    --recycle_early_stop_tolerance <float>: (default: 0.5) early stop tolerance for recycling.
    --use_dropout <bool>: (default: false) use dropout during inference.
    --stop_at_score <float>: (default: 100) stop at this confidence score.
    --disable_cluster_profile <bool>: (default: false) disable cluster profile.
    
Example in the environment of a virtual machine in the IFB-Biosphere cloud:
    nextflow main.nf --sequence examples/H1140.fasta --run test --database_dir ~/data/public/colabfold/ -profile docker -resume
'''

// -------------------------
// Paramètres par défaut (ligne de commande)
// -------------------------
params.help = false
params.h = false
params.sequence = null
params.run = null
params.database_dir = null
params.tool = 'ColabFold'
params.config_json = './massivefold_runs/ColabFold_params.json'
params.msas_precomputed = null
params.predictions_per_model = 1
params.batch_size = 5
params.pair_strategy = 'greedy'
params.num_recycle = 3
params.recycle_early_stop_tolerance = 0.5
params.use_dropout = false
params.stop_at_score = 100
params.disable_cluster_profile = false

// -------------------------
// Fonctions de configuration
// -------------------------

// Fonction pour charger toute la configuration
def loadFullConfig(jsonFile) {
    if (!file(jsonFile).exists()) {
        log.warn("Configuration file not found: ${jsonFile}. Using default/command-line parameters.")
        return [:]
    }
    def jsonSlurper = new groovy.json.JsonSlurper()
    return jsonSlurper.parse(file(jsonFile))
}

// Paramètres en commun
def extractCommonParams(config) {
    return [
        run_massivefold_plots      : config.massivefold?.run_massivefold_plots ?: '',
        data_dir                   : config.massivefold?.data_dir ?: '',
        jobfile_headers_dir        : config.massivefold?.jobfile_headers_dir ?: '',
        jobfile_templates_dir      : config.massivefold?.jobfile_templates_dir ?: '',
        scripts_dir                : config.massivefold?.scripts_dir ?: '',
        output_dir                 : config.massivefold?.output_dir ?: './output',
        logs_dir                   : config.massivefold?.logs_dir ?: './log',
        input_dir                  : config.massivefold?.input_dir ?: './input',
        models_to_use              : config.massivefold?.models_to_use ?: '',
        pkl_format                 : config.massivefold?.pkl_format ?: 'full',
        // Plots
        MF_plots_top_n_predictions : config.plots?.MF_plots_top_n_predictions ?: '10',
        MF_plots_chosen_plots      : config.plots?.MF_plots_chosen_plots ?: ''
    ]
}

// Spécifiques à chaque outil
def applyToolConfig(config, tool) {
    def toolParams = extractCommonParams(config)  // base commune

    switch(tool.toLowerCase()) {
        case 'colabfold':
            toolParams << [
                // CF_run params
                model_preset              : config.CF_run?.model_preset ?: 'multimer',
                pair_strategy             : config.CF_run?.pair_strategy ?: 'greedy',
                use_dropout               : config.CF_run?.use_dropout ?: 'false',
                num_recycle               : config.CF_run?.num_recycle ?: '20',
                recycle_early_stop_tolerance: config.CF_run?.recycle_early_stop_tolerance ?: '0.5',
                stop_at_score             : config.CF_run?.stop_at_score ?: '100',
                disable_cluster_profile   : config.CF_run?.disable_cluster_profile ?: 'false'
            ]
            break

        case 'afmassive':
            toolParams << [
                run_massivefold           : config.massivefold?.run_massivefold ?: '',
                uniref_database           : config.massivefold?.uniref_database ?: '',
                // AFM_run params
                model_preset              : config.AFM_run?.model_preset ?: 'multimer',
                max_recycles              : config.AFM_run?.max_recycles ?: '20',
                templates                 : config.AFM_run?.templates ?: 'true',
                dropout                   : config.AFM_run?.dropout ?: 'false',
                dropout_structure_module  : config.AFM_run?.dropout_structure_module ?: 'false',
                dropout_rates_filename    : config.AFM_run?.dropout_rates_filename ?: '',
                stop_recycling_below      : config.AFM_run?.stop_recycling_below ?: '0',
                max_template_date         : config.AFM_run?.max_template_date ?: '2024-08-31',
                min_score                 : config.AFM_run?.min_score ?: '0',
                max_score                 : config.AFM_run?.max_score ?: '1',
                db_preset                 : config.AFM_run?.db_preset ?: 'full_dbs',
                early_stop_tolerance      : config.AFM_run?.early_stop_tolerance ?: '0.5',
                bfd_max_hits              : config.AFM_run?.bfd_max_hits ?: '100000',
                mgnify_max_hits           : config.AFM_run?.mgnify_max_hits ?: '501',
                uniprot_max_hits          : config.AFM_run?.uniprot_max_hits ?: '50000',
                uniref_max_hits           : config.AFM_run?.uniref_max_hits ?: '10000'
            ]
            break

        case 'alphafold3':
            toolParams << [
                // Ajouter les paramètres spécifiques à AlphaFold3
                model_preset              : config.AF3_run?.model_preset ?: 'multimer',
                // Autres paramètres AF3...
            ]
            break

        default:
            log.warn "Unknown tool: ${tool}, using default parameters"
    }

    return toolParams
}

// -------------------------
// Chargement et fusion de la configuration
// -------------------------

// Charger la configuration JSON si elle existe
def fullConfig = loadFullConfig(params.config_json)

// Appliquer la configuration spécifique à l'outil
if (fullConfig) {
    def toolConfig = applyToolConfig(fullConfig, params.tool)
    
    // Fusionner avec les paramètres de ligne de commande (ligne de commande prioritaire)
    toolConfig.each { k, v -> 
        if (!params.containsKey(k) || params[k] == null) {
            params[k] = v
        }
    }
}

// Créer le fichier de configuration temporaire pour les processus
// Note: Ce fichier sera créé dans le workflow et passé aux processus

// -------------------------
// Workflow principal
// -------------------------

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

    // Normaliser le nom de l'outil
    def toolNormalized = params.tool.toLowerCase()

    // Debug information
    log.info """
    ==============================================
    MassiveFold Pipeline - ${params.tool.toUpperCase()}
    ==============================================
    Tool:                    ${params.tool}
    Sequence:                ${params.sequence}
    Run name:                ${params.run}
    Database directory:      ${params.database_dir}
    Pair strategy:           ${params.pair_strategy}
    Predictions per model:   ${params.predictions_per_model}
    Batch size:              ${params.batch_size}
    Num recycle:             ${params.num_recycle}
    Use dropout:             ${params.use_dropout}
    Stop at score:           ${params.stop_at_score}
    ==============================================
    """.stripIndent()

    // Treatment of the Multi-Sequence Alignment (MSA)
    def msa_results
    def seqFiles = Channel.fromPath(params.sequence)

    // Check the Tool used and apply the appropriate workflow
    if (toolNormalized == 'colabfold') {
        
        // Créer le fichier de configuration pour les processus
        def configToolContent = fullConfig ? groovy.json.JsonOutput.toJson(fullConfig) : '{}'
        configToolFile = Channel.of(configToolContent)
            .collectFile(name: 'config_tool.json', newLine: false)
        
        fastafile_unified = Unifier_colabfold(seqFiles)
        
        if (params.msas_precomputed) {
            // Use precomputed alignments
            log.info("Using precomputed MSAs: ${params.msas_precomputed}")
            msa_results = Channel
                .fromPath("${params.msas_precomputed}/*")
                .map { precomputed_msa -> tuple(precomputed_msa.baseName, precomputed_msa) }
        } else {
            // Run alignment process
            log.info("Running MSA alignment process")
            msa_results = Alignement_with_colabfold(
                fastafile_unified, 
                params.database_dir, 
                params.pair_strategy
            )
        }

        // Create batches
        batch_json = Create_batchs(
            fastafile_unified, 
            params.run, 
            params.predictions_per_model, 
            params.batch_size, 
            configToolFile, 
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
        out = Standardize_output_colabfold(
            fastafile_unified, 
            batch_json, 
            res_prediction.collect()
        )

    } 
    else if (toolNormalized == 'afmassive') {
        log.error("AFMassive tool is not yet implemented")
        exit 1
    }
    else if (toolNormalized == 'alphafold3') {
        fastafile_unified = Unifier_AlphaFold3(seqFiles)
        
        if (params.msas_precomputed) {
            log.info("Using precomputed MSAs: ${params.msas_precomputed}")
            msa_results = Channel
                .fromPath("${params.msas_precomputed}/*")
                .map { precomputed_msa -> tuple(precomputed_msa.baseName, precomputed_msa) }
        } else {
            msa_results = Alignement_with_colabfold(
                fastafile_unified, 
                params.database_dir, 
                params.pair_strategy
            )
        }
        
        log.error("AlphaFold3 tool is not yet fully implemented")
        exit 1
    } 
    else {
        log.error("Unsupported tool: ${params.tool}")
        log.error("Supported tools: ColabFold, AFMassive, AlphaFold3")
        exit 1
    }
}

// -------------------------
// Processus
// -------------------------

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

process Unifier_AlphaFold3 {
    label 'python_treatment'

    input:
    path sequence

    output:
    path "input/alphafold3_json_requests/${sequence.baseName}.fasta"

    script:
    """
    mkdir -p input/alphafold3_json_requests
    cp ${sequence} input/.
    unifier.py --conversion input --to_convert input/${sequence.baseName}.fasta --tool AlphaFold3 --json_params AlphaFold3_params.json
    """

    stub:
    """
    mkdir -p input/alphafold3_json_requests
    cp ${sequence} input/alphafold3_json_requests/${sequence.baseName}.fasta
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

    # Check if data directory exists
    if [ ! -d "${data_dir}" ]; then
        echo "ERROR: Database directory does not exist: ${data_dir}"
        exit 1
    fi
    
    # Set pairing strategy
    if [[ "${pair_strategy}" == "greedy" ]]; then
        pairing_strategy=0
    elif [[ "${pair_strategy}" == "complete" ]]; then
        pairing_strategy=1
    else
        echo "ValueError: --pair_strategy '${pair_strategy}' is not valid. Use 'greedy' or 'complete'"
        exit 1
    fi

    echo "Using pairing strategy: \${pairing_strategy}"
    
    # Run ColabFold search
    colabfold_search ${seqFile} ${data_dir} ${seqFile.baseName}_msa --pairing_strategy \${pairing_strategy}
    
    echo "=== ColabFold alignment completed ==="
    """

    stub:
    """
    echo "=== Running in stub mode ==="
    msa_dir="${seqFile.baseName}_msa"
    mkdir -p "\${msa_dir}"

    grep '^>' ${seqFile} | while read -r line; do
        header=\$(echo "\${line}" | sed 's/^>//')
        sanitized_header=\$(echo "\${header}" | tr -c '[:alnum:]_' '_')
        output_file="\${msa_dir}/${seqFile.baseName}_\${sanitized_header}.a3m"
        echo "Stub MSA content for \${line}" > "\${output_file}"
        echo "Created stub file: \${output_file}"
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
    # Detect GPU availability
    if command -v nvidia-smi &> /dev/null && nvidia-smi -L | grep -q "GPU"; then
        echo "✅ GPU detected, using CUDA"
        export CUDA_VISIBLE_DEVICES=0
    else
        echo "⚠️ No GPU detected, falling back to CPU"
        export CUDA_VISIBLE_DEVICES=""
        export JAX_PLATFORMS=cpu
    fi

    echo "=== Starting ColabFold inference ==="
    echo "Batch ID: ${id_batch}"
    echo "Sequence: ${sequence_name}"
    echo "Batch range: ${batch_start} to ${batch_end}"
    echo "Model: ${batch_model}"
    echo "MSA folder: ${msaFolder}"
    
    # Generate random seed
    random_seed=\$(shuf -i 0-1000000 -n 1)
    echo "Random seed: \${random_seed}"
    
    # Parse model information
    num_models=\$(echo ${batch_model} | cut -d "_" -f 2)
    echo "Model number: \${num_models}"
    
    num_version=\$(echo ${batch_model} | cut -d "v" -f 2)
    version_type="alphafold2_multimer_v\${num_version}"
    model_type=\${version_type}
    echo "Model type: \${model_type}"
    
    # Calculate number of seeds
    num_seeds=\$((${batch_end} - ${batch_start} + 1))
    echo "Number of predictions: \${num_seeds}"
    
    # Check MSA folder
    if [ ! -d "${msaFolder}" ]; then
        echo "ERROR: MSA folder does not exist: ${msaFolder}"
        exit 1
    fi
    
    echo "MSA folder contents:"
    ls -la ${msaFolder}
    
    # Prepare additional flags
    additional_flags=""
    if [ "${use_dropout}" = "true" ]; then
        additional_flags="\${additional_flags} --use-dropout"
    fi
    if [ "${disable_cluster_profile}" = "true" ]; then
        additional_flags="\${additional_flags} --disable-cluster-profile"
    fi
    
    # Run ColabFold batch prediction
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
    
    echo "=== ColabFold inference completed ==="
    """

    stub:
    """
    echo "=== Running inference in stub mode ==="
    output_dir="res_${sequence_name}_${run_name}_${id_batch}"
    mkdir -p "\${output_dir}"
    
    echo "Stub prediction output" > "\${output_dir}/prediction_stub.txt"
    echo "Model: ${batch_model}" > "\${output_dir}/model_info.txt"
    echo "Batch range: ${batch_start} to ${batch_end}" > "\${output_dir}/batch_info.txt"
    
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
    echo "Sequence file: ${seqFile}"
    echo "Batch JSON: ${batchjson}"
    
    mkdir -p organized_outputs
    
    # Run unifier
    unifier.py --conversion output --to_convert . --batches_file ${batchjson} --tool ColabFold
    
    # Organize outputs
    organize_outputs.py --batches_path .
    
    # Move organized outputs
    if [ -d "organized" ]; then
        mv organized/* organized_outputs/
    else
        cp -r res_* organized_outputs/ 2>/dev/null || true
    fi
    
    echo "=== Output standardization completed ==="
    """

    stub:
    """
    echo "=== Standardizing outputs in stub mode ==="
    mkdir -p organized_outputs
    echo "Stub standardized output" > organized_outputs/standardized_stub.txt
    echo "=== Stub standardization completed ==="
    """
}