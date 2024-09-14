# RNA-seq Pipeline for HPC
This repository contains a comprehensive RNA-seq analysis pipeline, designed to be executed on a High-Performance Computing (HPC) environment using SLURM for job scheduling. The pipeline processes RNA-seq data from raw reads to downstream analysis, including alternative splicing analysis using RMATS. 

## Features
Download RNA-seq data from SRA
Convert SRA to FASTQ format
Quality control and trimming
Alignment to reference genome (STAR or your preferred aligner)
BAM file sorting and indexing
Differential alternative splicing analysis using RMATS
Fully containerized (using Singularity or Docker)
Supports parallelized and efficient job execution on HPC environments