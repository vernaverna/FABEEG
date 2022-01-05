#!/bin/bash

# Make sure to request only the resources you really need to avoid cueing
#SBATCH --time=00-30
#SBATCH --mem-per-cpu=10000
#SBATCH -n 1
#SBATCH --output=./slurm_log/slurm-%A_%a.out

# Run for each subject (in this case, 9 times)
#SBATCH --array=1-9
INPUT_PATH=/m/nbe/scratch/rubberboot/BIDS/

# define n
n=$SLURM_ARRAY_TASK_ID                  

# get n:th line (1-indexed) of the file
SUB_ID=`sed "${n}q;d" aging_subjects.txt`    

# Load the Python environment. 
ml anaconda3
SUBJECTS_DIR=/m/nbe/scratch/rubberboot/BIDS/mri/
CODE_DIR=/m/nbe/scratch/rubberboot/code/

# Run the analysis
srun python ${CODE_DIR}/epochs.py $SUB_ID

