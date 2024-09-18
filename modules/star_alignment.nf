process star_alignment {
    publishDir "bam_files", mode: 'copy'
    publishDir "count_files", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path 'bam_files/*.bam', emit: bam_ch
    path 'count_files/*.out.tab'

    cpus params.cpus
    memory params.memory
    queue params.queue

    script:
    """
    module purge
    module load ${params.STAR}

    mkdir -p bam_files count_files

    STAR --genomeDir "${params.genomedir}" \
         --readFilesIn ${reads} \
         --readFilesCommand zcat \
         --quantMode GeneCounts \
         --outSAMtype BAM Unsorted \
         --sjdbGTFfile "${params.gtf_path}" \
         --outFileNamePrefix "${sample_id}_"

    mv *Aligned.out.bam bam_files/${sample_id}.bam
    mv *ReadsPerGene.out.tab count_files/${sample_id}.out.tab
    """
}
