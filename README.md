# RNA-seq Pipeline for HPC

This repository contains a comprehensive RNA-seq analysis pipeline tailored for High-Performance Computing (HPC) environments, utilizing SLURM for job scheduling. The pipeline processes RNA-seq data from raw reads to downstream analyses, including alternative splicing using RMATS.

## Features

- **Download RNA-seq data** from the SRA database.
- **Convert SRA files to FASTQ format**.
- **Quality control and trimming** of raw reads.
- **Alignment** to a reference genome (using STAR or other preferred aligners).
- **BAM file sorting and indexing**.
- **Differential alternative splicing analysis** using RMATS.
- **Containerized execution** using Singularity or Docker.
- **Parallelized and efficient job execution** designed for HPC environments.

## SLURM Job Support

Each process in the pipeline can be run individually as a SLURM job, allowing for flexibility in the execution of RNA-seq tasks. Individual SLURM sbatch scripts are provided for each step, giving users control over the workflow.

## Nextflow Integration

In addition to individual scripts, the repository includes a modular **Nextflow script** that enables users to select which processes to run. The Nextflow pipeline is designed to be executed as a SLURM job, providing:
- A streamlined and automated workflow.
- Modular execution, allowing users to choose specific processes.
- The ability to harness the power of parallel processing with Nextflow.

## Usage

### Running Individual Scripts

To run individual processes as SLURM jobs, follow these steps:

1. **Prepare SLURM sbatch Scripts**: Review and adjust the provided SLURM sbatch scripts for your HPC environment, ensuring the parameters match your computational resources and requirements. All slurm scripts are in slurm_scripts directory.

2. **Submit SLURM Jobs**: Use the `sbatch` command to submit each process to the SLURM queue. For example:
   ```bash
   sbatch sra_download.sh
   sbatch star_alignment.sh
   sbatch bam_sort.sh
   sbatch rmats.sh

### Running the Nextflow Pipeline

To run the entire pipeline using Nextflow with SLURM, follow these steps:

1. **Prepare Environment**: Ensure you have Nextflow installed and configured, and that your HPC environment supports Nextflow and SLURM. You may also need to load or activate the Nextflow environment. For example:

    ```bash
    source activate nextflow_env

2. **Edit SLURM Job Script**: Customize the rnaseq_pipeline_slurm.sh script with your specific parameters and resource settings.

3. **Submit Nextflow Job**: Submit the job to the SLURM queue using sbatch:

    ```bash
    sbatch rnaseq_pipeline_slurm.sh

## Parameters

- `-entry [optional] workflow_name`: 
  Select the specific workflow to execute (e.g., `do_full`, `do_fastq2bam`, `do_rmats`, `do_rmatsv2`).

- `-cartpath "path/to/cartfile"`: 
  Path to the SRA cartfile for data download.

- `-ngcpath "path/to/ngcfile"`: 
  Path to the NGC file for SRA authentication.

- `-genomedir "path/to/ref/genome/directory"`: 
  Path to the reference genome directory for alignment.

- `-gtfpath "path/to/annotated/gtf"`: 
  Path to the GTF file for gene annotation.

- `-rmats_b1 "path/to/rmats/batch1/files/list/txt"`: 
  Path to the list of batch 1 files for RMATS analysis.

- `-rmats_b2 "path/to/rmats/batch2/files/list/txt"`: 
  Path to the list of batch 2 files for RMATS analysis.

- `-singularity_image "path/to/singularity/image/for/rmats"`: 
  Path to the Singularity image for RMATS analysis.

- `-modulepath "path/to/modulefiles"`: 
  Path to modulefiles directory.

- `-sratoolkit "sratoolkit/[module version OPTIONAL]"`: 
  SRA Toolkit module.

- `-STAR "STAR/[module version OPTIONAL]"`: 
  STAR aligner module.

- `-samtools "samtools/[module version OPTIONAL]"`: 
  Samtools module.

- `-singularity "singularity/[module version OPTIONAL]"`: 
  Singularity module.

- `-intel "intel/[module version OPTIONAL]"`: 
  Intel module.

- `-pgi "pgi/[module version OPTIONAL]"`: 
  PGI module.

- `-rmats "rmats/[module version OPTIONAL]"`: 
  RMATS module.

- `-workingdir $(PWD)`: 
  Working directory for intermediate and output files.

- `-cpus ${SLURM_CPUS_PER_TASK}`: 
  Number of CPUs allocated for the job.

- `-memory ${SLURM_MEM}`: 
  Amount of memory allocated for the job.

- `-queue 'your_queue'`: 
  SLURM queue to use for job submission.



# Contribution

Feel free to contribute to this project by opening issues or submitting pull requests. Contributions are welcome to improve functionality, add features, or fix bugs.


# License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

