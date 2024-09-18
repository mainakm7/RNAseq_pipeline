#!/usr/bin/env bash

#SBATCH --partition=PARTITION_NAME    # Replace with your partition name
#SBATCH --requeue                     # Return job to the queue if preempted
#SBATCH --job-name=RNAseq_pipeline    # Assign a name to your job
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks=1                    # Total # of tasks across all nodes
#SBATCH --cpus-per-task=64            # Number of cores per task
#SBATCH --mem=200G                    # Amount of RAM required
#SBATCH --time=120:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out      # STDOUT output file
#SBATCH --error=slurm.%N.%j.err       # STDERR output file
#SBATCH --export=ALL                  # Export current environment to job
#SBATCH --mail-type=END,FAIL          # Send email on end or failure
#SBATCH --mail-user=EMAIL_ID          # Replace with your email address


source activate nextflow_env


nextflow run rnaseq_main.nf \
  -entry [optional] workflow_name \    # Specify workflow entry if needed
  -cartpath "path/to/cartfile" \
  -ngcpath "path/to/ngcfile" \
  -genomedir "path/to/ref/genome/directory" \
  -gtf_path "path/to/annotated/gtf" \
  -rmats_b1 "path/to/rmats/batch1/files/list/txt" \
  -rmats_b2 "path/to/rmats/batch2/files/list/txt" \
  -singularity_image "path/to/singularity/image/for/rmats" \
  -modulepath "path/to/modulefiles" \
  -sratoolkit "sratoolkit/[module version OPTIONAL]" \
  -STAR "STAR/[module version OPTIONAL]" \
  -samtools "samtools/[module version OPTIONAL]" \
  -singularity "singularity/[module version OPTIONAL]" \
  -intel "intel/[module version OPTIONAL]" \
  -pgi "pgi/[module version OPTIONAL]" \
  -rmats "rmats/[module version OPTIONAL]" \
  -workingdir $(PWD) \
  -cpus ${SLURM_CPUS_PER_TASK} \
  -memory ${SLURM_MEM} \
  -queue 'your_queue'   # Replace with your queue name
