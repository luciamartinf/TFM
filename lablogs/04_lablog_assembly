#!/usr/bin/env bash
#SBATCH --job-name metaspades       # Name for your job
#SBATCH --partition=bigmem,long     # total run time limit in HH:MM:SS
#SBATCH --time=99:00:00             # total run time limit in HH:MM:SS
#SBATCH --output=logs/metaspades-%a.%A.out # STDOUT file
#SBATCH --error=logs/metaspades-%a.%A.err  # STDERR file
#SBATCH --mail-type=end             # send email when job ends
#SBATCH --mail-user=luciamartinfernandez99@gmail.com # this is the email you wish to be notified at
#SBATCH --mem=150GB                 # Reserve X GB RAM for the job
#SBATCH -c 15                       # request for cores to be in the same node 
#SBATCH --array=0-116               # Array range and number of simultanous jobs


# Remember to create logs directory in the terminal
# mkdir logs 

# Load required modules
ml SPAdes/3.15.2-GCC-8.2.0-2.31.1
module load MEGAHIT/1.2.8-GCCcore-8.2.0 

# Create variables of input and output data
RUN_PATH_I=/home/lmartin/SRV_RyC/02-trimmomatic
RUN_PATH_O=/home/lmartin/SRV_RyC/05-assembly
names=($(cat /home/lmartin/SRV_RyC/00-raw_reads/samples_id.txt))

# Run assembly of each sample
# cat /home/sbodi/SRV_024/00-raw_reads/samples_id.txt | xargs -I % echo "mkdir -p %" > samples.sh

## METASPADES

/beegfs/easybuild/CentOS/7.6.1810/Skylake/software/SPAdes/3.15.2-GCC-8.2.0-2.31.1/bin/metaspades.py -1 $RUN_PATH_I/${names[${SLURM_ARRAY_TASK_ID}]}_R1.clean_qc_pair.fastq -2 $RUN_PATH_I/${names[${SLURM_ARRAY_TASK_ID}]}_R2.clean_qc_pair.fastq -k auto -t 8 -m 1000 -o $RUN_PATH_O/${names[${SLURM_ARRAY_TASK_ID}]}/

## MEGAHIT

# megahit -1 $RUN_PATH_I/${names[${SLURM_ARRAY_TASK_ID}]}_R1.clean_qc_pair.fastq -2 $RUN_PATH_I/${names[${SLURM_ARRAY_TASK_ID}]}_R2.clean_qc_pair.fastq --presets meta-large -o $RUN_PATH_O/${names[${SLURM_ARRAY_TASK_ID}]}/ -t 30 -m 1

# megahit -1 $RUN_PATH_I/R17W8_host_removed_R1.fastq.gz -2 $RUN_PATH_I/R17W8_host_removed_R2.fastq.gz --presets meta-large  -o $RUN_PATH_O/Megahit/R17W8





