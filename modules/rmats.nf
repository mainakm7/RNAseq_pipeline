process rmats {
    publishDir "${params.workingdir}/rmats_files", mode: 'copy'

    input:
    path batch1_file
    path batch2_file

    output:
    path "${params.workingdir}/rmats_files/*.txt"

    cpus params.cpus
    memory params.memory
    queue params.queue

    script:
    """
    module purge
    module load ${params.STAR}
    module load ${params.singularity}
    module load ${params.intel}
    module load ${params.pgi}

    mkdir -p ${params.workingdir}/rmats_files ${params.workingdir}/rmats_files/tmp

    singularity exec \
        -B ${workDir}:/data \
        ${params.singularity_image} \
        python -u /rmats-turbo/rmats.py \
            --b1 ${batch1_file} \
            --b2 ${batch2_file} \
            --gtf ${params.gtf_path} \
            -t paired \
            --readLength 100 \
            --od "${params.workingdir}/rmats_files" \
            --variable-read-length \
            --nthread ${task.cpus} \
            --tmp "${params.workingdir}/rmats_files/tmp" \
            --allow-clipping
    """
}