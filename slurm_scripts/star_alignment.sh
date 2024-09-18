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
module load STAR

ls *_1.fastq.gz > list.txt
sed -i 's/_1.fastq.gz//g' list.txt


GENOME_REF_DIR="/path/to/genomeDir"
GTF_PATH="/path/to/annotation.gtf"


while read -r line; do
    STAR --genomeDir "${GENOME_REF_DIR}" \
         --runThreadN ${SLURM_CPUS_PER_TASK} \
         --readFilesIn "${line}_1.fastq.gz" "${line}_2.fastq.gz" \
         --readFilesCommand zcat \
         --quantMode GeneCounts \
         --outSAMtype BAM Unsorted \
         --sjdbGTFfile "${GTF_PATH}" \
         --outFileNamePrefix "${line}"
done < list.txt
