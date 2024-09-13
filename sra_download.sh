#!/usr/bin/env bash

#SBATCH --partition=PARTITION_NAME    # Partition (job queue)
#SBATCH --requeue                     # Return job to the queue if preempted
#SBATCH --job-name=SRA_DOWNLOAD       # Assign a short name to your job
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks=1                    # Total # of tasks across all nodes
#SBATCH --cpus-per-task=64            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=200G                    # Real memory (RAM) required (MB)
#SBATCH --time=120:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out      # STDOUT output file
#SBATCH --error=slurm.%N.%j.err       # STDERR output file
#SBATCH --export=ALL                  # Export current environment to job
#SBATCH --mail-type=END,FAIL          # Send email on end or failure
#SBATCH --mail-user=EMAIL_ID          # Where to send email


module purge
module use /path/to/modulefiles/
module load sratoolkit/2.10.8

SRA_PROJECT_KEY="/path/to/project_key.ngc"
SRA_FILE_CART="/path/to/file_cart.krt"

srun prefetch -X 200G --ngc ${SRA_PROJECT_KEY} --cart ${SRA_FILE_CART}
