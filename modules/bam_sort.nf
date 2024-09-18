process bam_sort {
    publishDir "${params.workingdir}/sorted_bam_files", mode: 'copy'

    input:
    path bam_file

    output:
    path "${params.workingdir}/sorted_bam_files/*.sorted.bam", emit: sorted_bam_ch

    cpus params.cpus
    memory params.memory
    queue params.queue

    script:
    """
    module purge
    module load ${params.samtools}

    mkdir -p ${params.workingdir}/sorted_bam_files

    samtools sort -n -@ ${task.cpus} -m 2G ${bam_file} -o ${params.workingdir}/sorted_bam_files/${bam_file.baseName}.sorted.bam
    """
}