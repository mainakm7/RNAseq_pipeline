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
Each process in the pipeline can be run individually as a SLURM job, allowing for flexibility in the execution of RNA-seq tasks. These individual scripts are designed for easy submission to SLURM, giving users full control over the workflow.

## Nextflow Integration
In addition to individual scripts, the repository includes a modular **Nextflow script** that enables users to select which processes to run. The Nextflow pipeline is also designed to be executed as a SLURM job, providing:
- A streamlined and automated workflow.
- Modular execution, allowing users to choose specific processes.
- The ability to harness the power of parallel processing with Nextflow.

## Usage
1. **Run Individual Scripts**: Submit each process independently as a SLURM job using the provided SLURM sbatch scripts.
2. **Run Nextflow Pipeline**: Use the Nextflow script to run selected processes, with full support for SLURM job submission and resource management.
