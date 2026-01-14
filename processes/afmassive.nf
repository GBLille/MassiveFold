process Alignement_with_AFMassive {
    tag "$seqFile.baseName"
    publishDir "result/${seqFile.baseName}/alignment"
    label 'afmassive'
    containerOptions "-v ${data_dir}:${data_dir}:rw"

    input:
    path(seqFile)
    val(data_dir)
    path(params_weight)
    val(output_dir)

    output:
    tuple val(seqFile.baseName), path("${output_dir}/${seqFile.baseName}/msas")

    script:
    """
    echo "=== Starting AFMassive alignment ==="
    echo "Sequence file: ${seqFile}"
    echo "Output directory: ${output_dir}"
    echo "Database directory: ${data_dir}"
    echo "Weight directory: ${params_weight}"

    if [ ! -d "${data_dir}" ]; then
        echo "ERROR: Database directory does not exist: ${data_dir}"
        exit 1
    fi

    if [ ! -d "${params_weight}" ]; then
        echo "ERROR: Database directory does not exist: ${params_weight}"
        exit 1
    fi

    run_AFmassive.py \
        --alignments_only=true \
        --fasta_paths=${seqFile} \
        --output_dir=${output_dir} \
        --data_dir=${params_weight} \
        --use_gpu_relax=false \
        --max_template_date=2024-08-31 \
        --db_preset=reduced_dbs \
        --model_preset=multimer \
        --small_bfd_database_path=${data_dir}/small_bfd/bfd-first_non_consensus_sequences.fasta \
        --uniref90_database_path=${data_dir}/uniref90/uniref90.fasta \
        --mgnify_database_path=${data_dir}/mgnify/mgy_clusters_2022_05.fa \
        --template_mmcif_dir=${data_dir}/pdb_mmcif/mmcif_files \
        --pdb_seqres_database_path=${data_dir}/pdb_seqres/pdb_seqres.txt \
        --uniprot_database_path=${data_dir}/uniprot/uniprot.fasta \
        --obsolete_pdbs_path=${data_dir}/pdb_mmcif/obsolete.dat



    echo "=== AFMassive alignment completed ==="
    """
}
process Unifier_AFMassive {
    label 'python_treatment'

    input:
    path sequence

    output:
    path "input/*.fasta"

    script:
    def seqName = sequence.getBaseName()
    """
    mkdir -p input
    cp ${sequence} input/${seqName}.fasta
    """
}


process Run_inference_AFMassive {
    tag "$sequence_name | batch#$id_batch"
    label 'afmassive'
    maxForks 2

    input:
    tuple val(id_batch), val(batch_start), val(batch_end), val(batch_model), val(sequence_name), path(msaFolder), val(run_name), path(seqFile), val(data_dir), val(params_weight), val(output_dir)

    output:
    path "batch_${id_batch}/*"

    script:
    """
    # Detect GPU availability
    if command -v nvidia-smi &> /dev/null && nvidia-smi -L | grep -q "GPU"; then
        echo "GPU detected, using CUDA"
        export CUDA_VISIBLE_DEVICES=0
    else
        echo "No GPU detected, falling back to CPU"
        export CUDA_VISIBLE_DEVICES=""
    fi

    echo "=== Starting AFMassive inference ==="
    echo "Batch ID: ${id_batch}"
    echo "Sequence: ${sequence_name}"
    echo "Batch range: ${batch_start} to ${batch_end}"
    echo "Model: ${batch_model}"
    echo "MSA folder: ${msaFolder}"

    # Check MSA folder
    if [ ! -d "${msaFolder}" ]; then
        echo "ERROR: MSA folder does not exist: ${msaFolder}"
        exit 1
    fi

    # Create output directory
    mkdir -p batch_${id_batch}/${sequence_name}

    mv ${msaFolder} batch_${id_batch}/${sequence_name}/msas

    # Run AFMassive inference
    run_AFmassive.py \
        --use_precomputed_msas=true \
        --fasta_paths=${seqFile} \
        --output_dir=batch_${id_batch} \
        --start_prediction=${batch_start} \
        --end_prediction=${batch_end} \
        --models_to_use=${batch_model} \
        --data_dir=${params_weight} \
        --use_gpu_relax=false \
        --max_template_date=2024-08-31 \
        --db_preset=reduced_dbs \
        --model_preset=multimer \
        --small_bfd_database_path=${data_dir}/small_bfd/bfd-first_non_consensus_sequences.fasta \
        --uniref90_database_path=${data_dir}/uniref90/uniref90.fasta \
        --mgnify_database_path=${data_dir}/mgnify/mgy_clusters_2022_05.fa \
        --template_mmcif_dir=${data_dir}/pdb_mmcif/mmcif_files \
        --pdb_seqres_database_path=${data_dir}/pdb_seqres/pdb_seqres.txt \
        --uniprot_database_path=${data_dir}/uniprot/uniprot.fasta \
        --obsolete_pdbs_path=${data_dir}/pdb_mmcif/obsolete.dat

    echo "=== AFMassive inference completed ==="
    """
}


process Organize_output_AFMassive {
    publishDir "result/prediction/${seqFile.baseName}/"
    label 'python_treatment'

    input:
    path seqFile
    path batchjson
    path afmassive_batch_files

    output:
    path "organized_outputs/*"

    script:
    """
    echo "=== Organizing AFMassive outputs ==="
    echo "Sequence file: ${seqFile}"
    echo "Batch JSON: ${batchjson}"
    
    mkdir -p organized_outputs

    # Organize outputs (AFMassive outputs are already in AF2 standard format)
    organize_outputs.py --batches_path .

    # Move organized outputs
    if [ -d "organized" ]; then
        mv organized/* organized_outputs/
    else
        cp -r batch_* organized_outputs/ 2>/dev/null || true
    fi

    echo "=== Output organization completed ==="
    """

    stub:
    """
    echo "=== Organizing outputs in stub mode ==="
    mkdir -p organized_outputs
    echo "Stub organized output" > organized_outputs/organized_stub.txt
    echo "=== Stub organization completed ===" 
    """
}