nextflow.enable.dsl=2

params {
    cartpath = ""
    ngcpath = ""
    genomedir = ""
    gtf_path = ""
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
    rmats = 'rmats'
    cpus = 64
    memory = "200GB"
    queue = 'your_queue'
}

process sra_download {
    publishDir "fastq_files", mode: 'copy'

    input:
    path ngcpath
    path cartpath

    output:
    path 'fastq_files/*.{fq,fq.gz,fastq,fastq.gz}'

    cpus params.cpus
    memory params.memory
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

process rmats {
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
    module load ${params.singularity}
    module load ${params.intel}
    module load ${params.pgi}

    mkdir -p rmats_files rmats_files/tmp

    singularity exec \
        -B ${workDir}:/data \
        ${params.singularity_image} \
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

workflow do_full {
    take:
    ngcpath
    cartpath

    main:
    sra_download(ngcpath, cartpath)
    fastqpair_ch = Channel.fromFilePairs('fastq_files/*_{1,2}.{fq,fq.gz,fastq,fastq.gz}', size: 2)
    star_alignment(fastqpair_ch)
    bam_sort(star_alignment.out.bam_ch)
}

workflow do_fastq2bam {
    take:
    fastqpair_ch = Channel.fromFilePairs('fastq_files/*_{1,2}.{fq,fq.gz,fastq,fastq.gz}', size: 2)

    main:
    star_alignment(fastqpair_ch)
    bam_sort(star_alignment.out.bam_ch)
}

workflow do_rmats {
    take:
    batch1_file
    batch2_file

    main:
    rmats(batch1_file, batch2_file)
}

workflow do_rmatsv2 {
    take:
    batch1_file
    batch2_file

    main:
    rmatsv2(batch1_file, batch2_file)
}

workflow {
    do_full(params.ngcpath, params.cartpath)
}
