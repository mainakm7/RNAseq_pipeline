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
    workingdir = "/path/to/working/directory"
    cpus = 64
    memory = "200GB"
    queue = 'your_queue'
}


include { sra_download } from './modules/sra_download.nf'
include { star_alignment } from './modules/star_alignment.nf'
include { bam_sort } from './modules/bam_sort.nf'
include { rmats } from './modules/rmats.nf'
include { rmatsv2 } from './modules/rmatsv2.nf'


workflow do_full {
    take:
    ngcpath
    cartpath

    main:
    sra_download(ngcpath, cartpath)
    fastqpair_ch = Channel.fromFilePairs("${params.workingdir}/fastq_files/*_{1,2}.{fq,fq.gz,fastq,fastq.gz}", size: 2)
    star_alignment(fastqpair_ch)
    bam_sort(star_alignment.out.bam_ch)
}

workflow do_fastq2bam {
    take:
    fastqpair_ch = Channel.fromFilePairs("${params.workingdir}/fastq_files/*_{1,2}.{fq,fq.gz,fastq,fastq.gz}", size: 2)

    main:
    star_alignment(fastqpair_ch)
    bam_sort(star_alignment.out.bam_ch)
}

workflow do_rmats {
    main:
    if (params.rmats_b1 && params.rmats_b2) {
        rmats(params.rmats_b1, params.rmats_b2)
    } else {
        error "Batch files for RMATS are missing!"
    }
}


workflow do_rmatsv2 {
    main:
    if (params.rmats_b1 && params.rmats_b2) {
        rmats(params.rmats_b1, params.rmats_b2)
    } else {
        error "Batch files for RMATS are missing!"
    }
}

workflow {
    do_full(params.ngcpath, params.cartpath)
}