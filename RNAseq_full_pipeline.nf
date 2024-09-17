params {
    cartpath = ""
    ngcpath = ""
    genomedir = ""
    gtf_path = ""
    sra_download = true
    star_alignment = true
    bam_sort = true
    rmats = true
    rmats_b1 = ""
    rmats_b2 = ""
    singularity_image = ""
    modulepath = ''
    sratoolkit = 'sratoolkit/2.10.8'
    STAR = 'STAR'
    samtools = 'samtools/1.3.1'
    singularity = 'singularity/3.1.0'
    intel = 'intel/17.0.4'
    pgi = 'pgi'
    cpus = 64
    memory = "200GB"
    queue = 'your_queue'
}

process sra_download {
    publishDir "fastq_files", type: 'directory'

    input:
    path ngcpath from params.ngcpath
    path cartpath from params.cartpath

    output:
    path 'fastq_files/*.{fq,fq.gz,fastq,fastq.gz}' into fastq_ch

    cpus ${params.cpus}
    memory ${params.memory}
    queue params.queue

    script:
    """
    module purge
    module use ${params.modulepath}
    module load ${params.sratoolkit}

    mkdir -p fastq_files

    srun prefetch -X 200G --ngc ${ngcpath} --cart ${cartpath}
    
    mv *.fq.gz *.fq *.fastq *.fastq.gz fastq_files/
    """
}

fastqpair_ch = Channel.fromFilePairs( 'fastq_files/*_{1,2}.{fq,fq.gz,fastq,fastq.gz}', size: 2 )

process star_alignment {
    publishDir "bam_files", type: 'directory'
    publishDir "count_files", type: 'directory'

    input:
    tuple path fastq_file1, fastq_file2 from fastqpair_ch

    output:
    path 'bam_files/*.bam' into bam_ch
    path 'count_files/*.out.tab' into countfiles_ch

    cpus ${params.cpus}
    memory ${params.memory}
    queue params.queue

    script:
    """
    module purge
    module load ${params.STAR}

    base_name=${fastq_file1.baseName.replace('_1', '')}
    GENOME_REF_DIR=${params.genomedir}
    GTF_PATH=${params.gtf_path}

    mkdir -p bam_files count_files

    STAR --genomeDir "${GENOME_REF_DIR}" \
         --readFilesIn "${fastq_file1}" "${fastq_file2}" \
         --readFilesCommand zcat \
         --quantMode GeneCounts \
         --outSAMtype BAM Unsorted \
         --sjdbGTFfile "${GTF_PATH}" \
         --outFileNamePrefix "${base_name}_"

    mv *.bam bam_files/
    mv *.out.tab count_files/
    """
}

process bam_sort {
    publishDir "sorted_bam_files", type: 'directory'

    input:
    path bam_file from bam_ch

    output:
    path 'sorted_bam_files/*.sorted.bam' into sorted_bam_ch

    cpus ${params.cpus}
    memory ${params.memory}
    queue params.queue

    script:
    """
    module purge
    module load ${params.samtools}

    mkdir -p sorted_bam_files

    samtools sort -n -@ ${task.cpus} -m 2G ${bam_file} -o sorted_bam_files/\${bam_file.baseName}.sorted.bam
    """
}

batch1_ch = Channel.value(file(params.rmats_b1).list())
batch2_ch = Channel.value(file(params.rmats_b2).list())

process rmats {
    publishDir "rmats_files", type: 'directory'

    input:
    val batch1_files from batch1_ch
    val batch2_files from batch2_ch

    output:
    path 'rmats_files/*.txt' into rmats_ch

    cpus ${params.cpus}
    memory ${params.memory}
    queue params.queue

    script:
    """
    module purge
    module load ${params.STAR}
    module load ${params.singularity}
    module load ${params.intel}
    module load ${params.pgi}

    mkdir -p rmats_files/tmp

    SINGULARITY_IMAGE=${params.singularity_image}

    singularity exec \
        -B ${workDir}:/data \
        ${SINGULARITY_IMAGE} \
        python -u /rmats-turbo/rmats.py \
            --b1 ${batch1_files.join(',')} \
            --b2 ${batch2_files.join(',')} \
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


workflow rna_seq {
    // Conditional logic for each process
    if (params.sra_download) {
        sra_download()
    } else {
        fastq_ch = Channel.fromPath("fastq_files/*.{fq,fq.gz,fastq,fastq.gz}")
    }

    if (params.star_alignment) {
        star_alignment()
    } else {
        bam_ch = Channel.fromPath("bam_files/*.bam")
    }

    if (params.bam_sort) {
        bam_sort()
    } else {
        sorted_bam_ch = Channel.fromPath("sorted_bam_files/*.sorted.bam")
    }

    if (params.rmats && params.rmats_b1 && params.rmats_b2) {
        rmats()
    }
}
