#!/bin/bash
#SBATCH --job-name=hello_serial
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --partition=develop

./hello_serial
