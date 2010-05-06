#!/bin/bash
#SBATCH --job-name=cell$cd$
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --partition=serial

./cell$cd$
