params.cartpath = ""
params.ngcpath = ""
params.genomedir = ""
params.gtf_path = ""
params.sra_download = true
params.star_alignment = true
params.rmats_b1 = ""
params.rmats_b2 = ""
params.cpus = 64
params.memory = "200GB"

process sra_download {

    publishDir "fastq_files", type: 'directory'

    input:
    path ngcpath from params.ngcpath
    path cartpath from params.cartpath

    output:
    path 'fastq_files/*.{fq,fq.gz,fastq,fastq.gz}' into fastq_ch

    cpus params.cpus
    memory params.memory

    script:
    """
    module purge
    module use /path/to/modulefiles/
    module load sratoolkit/2.10.8

    mkdir -p fastq_files

    srun prefetch -X 200G --ngc ${ngcpath} --cart ${cartpath}
    
    mv *.fq.gz *.fq *.fastq *.fastq.gz fastq_files/
    """
}

process star_alignment {

    publishDir "bam_files", type: 'directory'
    publishDir "count_files", type: 'directory'

    input:
    path fastq_file from fastq_ch

    output:
    path 'bam_files/*.bam' into bam_ch
    path 'count_files/*.out.tab' into countfiles_ch

    cpus params.cpus
    memory params.memory

    script:
    """
    module purge
    module load STAR

    # Extract the base name from the fastq file
    base_name=\$(basename ${fastq_file} _1.fastq.gz)

    GENOME_REF_DIR=${params.genomedir}
    GTF_PATH=${params.gtf_path}

    mkdir -p bam_files count_files

    STAR --genomeDir "${GENOME_REF_DIR}" \
         --readFilesIn "${base_name}_1.fastq.gz" "${base_name}_2.fastq.gz" \
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

    cpus params.cpus
    memory params.memory

    script:
    """
    module purge
    module load samtools/1.3.1

    mkdir -p sorted_bam_files

    ls *.bam > list.txt

    while read -r line; do
        samtools sort -n -@ ${task.cpus} -m 2G ${line} -o ${line}.sorted.bam
    done < list.txt

    mv *.sorted.bam sorted_bam_files/
    """
}

process rmats {

    publishDir "rmats_files", type: 'directory'

    input:
    path sorted_bam_file from sorted_bam_ch

    output:
    path 'rmats_files/*.txt' into rmats_ch

    cpus params.cpus
    memory params.memory

    script:
    """
    module purge
    module load STAR/2.7.5a
    module load singularity/3.1.0
    module load intel/17.0.4
    module load pgi

    BATCH1_PATH=${params.rmats_b1}
    BATCH2_PATH=${params.rmats_b2}

    # Create results directory if it doesn't exist
    mkdir -p rmats_files/tmp

    # Define Singularity image
    SINGULARITY_IMAGE="/path/to/your_image.sif"

    # Run RMATS on both batches in one Singularity command
    singularity exec \
        -B ${workDir}:/data \
        ${SINGULARITY_IMAGE} \
        python -u /rmats-turbo/rmats.py \
            --b1 ${BATCH1_PATH} \
            --b2 ${BATCH2_PATH} \
            --gtf /path/to/annotated.gtf \
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
    if (params.sra_download) {
        sra_download()
    }

    if (params.star_alignment) {
        star_alignment()
    }

    if (params.bam_sort) {
        bam_sort()
    }

    if (params.rmats) {
        rmats()
    }
}
