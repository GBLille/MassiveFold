nextflow.enable.dsl = 2
import groovy.json.JsonSlurper

include {
    Unifier_colabfold;
    Alignement_with_colabfold;
    Run_inference_colabfold;
    Standardize_output_colabfold
} from './processes/colabfold.nf'

include {
    Alignement_with_AFMassive;
    Unifier_AFMassive;
    Run_inference_AFMassive;
    Organize_output_AFMassive
} from './processes/afmassive.nf'

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
// Workflow principal
// -------------------------

workflow {
    // Display the help message if requested
    if (params.help || params.h) {
        log.info(helpMessage)
        exit 0
    }

    // Validate required parameters
    if (!params.sequence || !params.run || !params.database_dir || !params.tool) {
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
    ==============================================
    """.stripIndent()

    // Vérifier que le bloc existe bien dans params
    if (!params[toolNormalized]) {
        log.error "Aucun paramètre spécifique trouvé pour l'outil '${params.tool}'."
    } else {
        log.info "Paramètres spécifiques à l'outil '${params.tool.toUpperCase()}':"
        params[toolNormalized].each { key, value ->
            log.info "   ${key} = ${value}"
        }
    }
    
    // Treatment of the Multi-Sequence Alignment (MSA)
    def msa_results
    def seqFiles = Channel.fromPath(params.sequence)

    // COLABFOLD 

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

        // Créer le fichier JSON des batches
        batch_json = Create_batches_json(
            fastafile_unified, 
            params.run, 
            params.predictions_per_model, 
            params.batch_size, 
            params.tool
        )
        
        // Parser le JSON et créer les tuples pour les batches (sans le json_file)
        batchs = batch_json
            .map { json_file ->
                def slurper = new JsonSlurper()
                def batches_data = slurper.parse(json_file)
                return batches_data
            }
            .flatMap { batches_data ->
                batches_data.collect { batch_id, batch_info ->
                    tuple(
                        batch_id as Integer,
                        batch_info.start as Integer,
                        batch_info.end as Integer,
                        batch_info.model
                    )
                }
            }

        // Combiner les batches avec les paramètres constants, puis avec les MSA results
        batches_with_params = batchs
            .map { batch_id, start, end, model ->
                tuple(batch_id, start, end, model, params.run, params.num_recycle, params.recycle_early_stop_tolerance, params.use_dropout, params.stop_at_score, params.disable_cluster_profile)
            }
        
        batches_msa = batches_with_params.combine(msa_results)

        // Run inference
        res_prediction = Run_inference_colabfold(batches_msa)
        
        // Standardize output
        out = Standardize_output_colabfold(
            fastafile_unified, 
            batch_json, 
            res_prediction.collect()
        )

    } 

    // AFMASSIVE


    else if (toolNormalized == 'afmassive') {
        // Étape 1 : unifier les séquences
        fastafile_unified = Unifier_AFMassive(seqFiles)

        // Étape 2 : utiliser des MSAs pré-calculés si dispo, sinon lancer AFMassive pour alignement
        if (params.msas_precomputed) {
            log.info("Using precomputed MSAs: ${params.msas_precomputed}")
            msa_results = Channel
                .fromPath("${params.msas_precomputed}")
                .map { msa_dir ->
                    tuple(
                        msa_dir.toFile().parentFile.name,
                        msa_dir
                    )
                }

            msa_results.view { v -> "DEBUG final tuple: $v" }

        } else {
            log.info("Running AFMassive alignment process")
            msa_results = Alignement_with_AFMassive(
                fastafile_unified,
                params.database_dir,
                params.afmassive.params_weight,
                params.output_dir
            )
        }

        // Étape 3 : création du JSON des batches
        batch_json = Create_batches_json(
            fastafile_unified,
            params.run,
            params.predictions_per_model,
            params.batch_size,
            params.tool
        )
        
        // Parser le JSON et créer les tuples (sans le json_file)
        batchs = batch_json
            .map { json_file ->
                def slurper = new JsonSlurper()
                def batches_data = slurper.parse(json_file)
                return batches_data
            }
            .flatMap { batches_data ->
                batches_data.collect { batch_id, batch_info ->
                    tuple(
                        batch_id as Integer,
                        batch_info.start as Integer,
                        batch_info.end as Integer,
                        batch_info.model
                    )
                }
            }

        // Étape 4 : combiner les batchs et MSAs

        batches_msa = batchs
            .combine(msa_results)
            .combine(fastafile_unified)
            .map { batch_id, start, end, model, sequence_name, msaFolder, fasta_file ->
                tuple(
                    batch_id,                       // id_batch
                    start,                          // batch_start
                    end,                            // batch_end
                    model,                          // batch_model
                    sequence_name,                  // sequence_name
                    msaFolder,                      // msaFolder (path)
                    params.run,                     // run_name
                    fasta_file,                     // seqFile (path)
                    params.database_dir,            // data_dir
                    params.afmassive.params_weight, // params_weight
                    params.output_dir               // output_dir
                )
            }
        msa_results.view  { v -> "msa: $v" }
        fastafile_unified.view  { v -> "fasta: $v" }
        batches_msa.view { v -> "Batch ready for inference: $v" }

        // Étape 5 : lancer les inférences
        res_prediction = Run_inference_AFMassive( batches_msa )

        // // Étape 6 : organiser les résultats finaux
        // out = Organize_output_AFMassive(
        //     fastafile_unified,
        //     batch_json,
        //     res_prediction
        // )
    }
    else if (toolNormalized == 'alphafold3') {
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
// Process globaux
// -------------------------

/**
 * Processus unique qui crée le JSON des batches
 * Retourne le fichier JSON qui sera parsé par Nextflow
 */
process Create_batches_json {
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
    
    // Récupérer le model_preset depuis les paramètres
    def modelPreset = toolCode == 'colabfold' ? params.colabfold?.model_preset ?: 'multimer' :
                      toolCode == 'afmassive' ? params.afmassive?.model_preset ?: 'multimer' :
                      'multimer'
    
    // Récupérer les modèles à utiliser
    def modelsToUse = params.massivefold?.models_to_use ?: []
    def modelsFilter = modelsToUse ? modelsToUse.join(',') : ''
    
    """
    python3 << 'PYEOF'
import math
import json

# Configuration
predictions_per_model = ${predictions_per_model}
batch_size = ${batch_size}
tool = "${toolCode}"
model_preset = "${modelPreset}"
models_filter = "${modelsFilter}".split(',') if "${modelsFilter}" else []

# Définir les modèles disponibles selon le preset
if tool == "alphafold3":
    model_names = ["AlphaFold3"]
elif model_preset == 'multimer':
    model_names = [
        'model_1_multimer_v1', 'model_2_multimer_v1', 'model_3_multimer_v1', 
        'model_4_multimer_v1', 'model_5_multimer_v1',
        'model_1_multimer_v2', 'model_2_multimer_v2', 'model_3_multimer_v2',
        'model_4_multimer_v2', 'model_5_multimer_v2',
        'model_1_multimer_v3', 'model_2_multimer_v3', 'model_3_multimer_v3',
        'model_4_multimer_v3', 'model_5_multimer_v3'
    ]
elif model_preset == 'monomer_ptm':
    model_names = [
        'model_1_ptm', 'model_2_ptm', 'model_3_ptm', 'model_4_ptm', 'model_5_ptm'
    ]
else:
    model_names = []

# Filtrer selon models_filter si défini
if models_filter and models_filter[0]:
    model_names = [m for m in model_names if m in models_filter]

print(f"Running inference on models: {', '.join(model_names)}")
print(f"Running {predictions_per_model} predictions on each of the {len(model_names)} models")
print(f"Total prediction number: {predictions_per_model * len(model_names)}")

# Calculer les batches par modèle
opt_batch_nb = math.ceil(predictions_per_model / batch_size)
batch_sizes = []
remaining = predictions_per_model

for _ in range(opt_batch_nb):
    if remaining >= batch_size:
        batch_sizes.append(batch_size)
        remaining -= batch_size
    else:
        batch_sizes.append(remaining)

# Calculer les edges des batches
batch_edges = [0]
for size in batch_sizes:
    batch_edges.append(batch_edges[-1] + size)

# Créer le dictionnaire JSON avec tous les batches
batches_dict = {}
batch_id = 0

for model in model_names:
    for i in range(len(batch_sizes)):
        start = batch_edges[i]
        end = batch_edges[i + 1] - 1
        batches_dict[str(batch_id)] = {
            'start': start,
            'end': end,
            'model': model
        }
        batch_id += 1

# Écrire le fichier JSON
output_file = "${sequence.baseName}_${run}_batches.json"
with open(output_file, 'w') as f:
    json.dump(batches_dict, f, indent=4)

print(f"Created {batch_id} batches in {output_file}")
PYEOF
    """
    
    stub:
    """
    cat > ${sequence.baseName}_${run}_batches.json << 'EOF'
{
    "0": {
        "start": 0,
        "end": 4,
        "model": "model_1_multimer_v3"
    }
}
EOF
    """
}

// process Unifier_AlphaFold3 {
//     label 'python_treatment'

//     input:
//     path sequence

//     output:
//     path "input/alphafold3_json_requests/${sequence.baseName}.fasta"

//     script:
//     """
//     mkdir -p input/alphafold3_json_requests
//     cp ${sequence} input/.
//     unifier.py --conversion input --to_convert input/${sequence.baseName}.fasta --tool AlphaFold3 --json_params AlphaFold3_params.json
//     """

//     stub:
//     """
//     mkdir -p input/alphafold3_json_requests
//     cp ${sequence} input/alphafold3_json_requests/${sequence.baseName}.fasta
//     """
// }