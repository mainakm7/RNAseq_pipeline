process sra_download {
    publishDir "${params.workingdir}/fastq_files", mode: 'copy'

    input:
    path ngcpath
    path cartpath

    output:
    path "${params.workingdir}/fastq_files/*.{fq,fq.gz,fastq,fastq.gz}"

    cpus params.cpus
    memory params.memory
    queue params.queue

    script:
    """
    module purge
    module use ${params.modulepath}
    module load ${params.sratoolkit}

    mkdir -p ${params.workingdir}/fastq_files

    srun prefetch -X 200G --ngc ${ngcpath} --cart ${cartpath}
    
    mv *.fq.gz *.fq *.fastq *.fastq.gz ${params.workingdir}/fastq_files/
    """
}