#!/bin/bash
#SBATCH --job-name=generic_quick
#SBATCH --time=04:00:00                    
#SBATCH --partition=4hrs
#SBATCH --ntasks=1
#SBATCH --mem=30G

srun python find_largest_burst.py\
    --input "$1"\
