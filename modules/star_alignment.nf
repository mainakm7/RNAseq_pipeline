process star_alignment {
    publishDir "${params.workingdir}/bam_files", mode: 'copy'
    publishDir "${params.workingdir}/count_files", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "${params.workingdir}/bam_files/*.bam", emit: bam_ch
    path "${params.workingdir}/count_files/*.out.tab"

    cpus params.cpus
    memory params.memory
    queue params.queue

    script:
    """
    module purge
    module load ${params.STAR}

    mkdir -p ${params.workingdir}/bam_files ${params.workingdir}/count_files

    STAR --genomeDir "${params.genomedir}" \
         --readFilesIn ${reads} \
         --readFilesCommand zcat \
         --quantMode GeneCounts \
         --outSAMtype BAM Unsorted \
         --sjdbGTFfile "${params.gtf_path}" \
         --outFileNamePrefix "${sample_id}_"

    mv *Aligned.out.bam ${params.workingdir}/bam_files/${sample_id}.bam
    mv *ReadsPerGene.out.tab ${params.workingdir}/count_files/${sample_id}.out.tab
    """
}