process rmatsv2 {
    publishDir "rmats_files", mode: 'copy'

    input:
    path batch1_file
    path batch2_file

    output:
    path 'rmats_files/*.txt'

    cpus params.cpus
    memory params.memory
    queue params.queue

    script:
    """
    module purge
    module load ${params.STAR}
    module load ${params.intel}
    module load ${params.pgi}

    module use ${params.modulepath}
    module load ${params.rmats}

    mkdir -p rmats_files rmats_files/tmp

    python -u /rmats-turbo/rmats.py \
        --b1 ${batch1_file} \
        --b2 ${batch2_file} \
        --gtf ${params.gtf_path} \
        -t paired \
        --readLength 100 \
        --od "rmats_files" \
        --variable-read-length \
        --nthread ${task.cpus} \
        --tmp "rmats_files/tmp" \
        --allow-clipping
    """
}
