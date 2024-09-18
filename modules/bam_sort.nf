process bam_sort {
    publishDir "sorted_bam_files", mode: 'copy'

    input:
    path bam_file

    output:
    path 'sorted_bam_files/*.sorted.bam', emit: sorted_bam_ch

    cpus params.cpus
    memory params.memory
    queue params.queue

    script:
    """
    module purge
    module load ${params.samtools}

    mkdir -p sorted_bam_files

    samtools sort -n -@ ${task.cpus} -m 2G ${bam_file} -o sorted_bam_files/${bam_file.baseName}.sorted.bam
    """
}
