#!/usr/bin/env bash

#SBATCH --partition=PARTITION_NAME    # Partition (job queue)
#SBATCH --requeue                     # Return job to the queue if preempted
#SBATCH --job-name=STARMAP            # Assign a short name to your job
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

# source activate conda_env

module load STAR
module load singularity
module load intel/17.0.4
module load pgi

WORKDIR="/path/to/working/directory"
BATCH1_PATH="/path/to/batch1_filenames.txt"
BATCH2_PATH="/path/to/batch2_filenames.txt"
RESULTS_DIR="${WORKDIR}/results"

# Create results directory if it doesn't exist
mkdir -p ${RESULTS_DIR}

# Define Singularity image
SINGULARITY_IMAGE="/path/to/your_image.sif"

# Run RMATS on both batches in one Singularity command
singularity exec \
    -B ${WORKDIR}:/data \
    ${SINGULARITY_IMAGE} \
    python -u /rmats-turbo/rmats.py \
        --b1 ${BATCH1_PATH} \
        --b2 ${BATCH2_PATH} \
        --gtf /path/to/annotated.gtf \
        -t paired \
        --readLength 100 \
        --od ${RESULTS_DIR} \
        --variable-read-length \
        --nthread ${SLURM_CPUS_PER_TASK} \
        --tmp ${RESULTS_DIR}/tmp \
        --allow-clipping

input > output
sacct --format MaxRSS,Elapsed -j $SLURM_JOBID