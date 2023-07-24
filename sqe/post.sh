#!/bin/bash -l
#SBATCH -C cpu
#SBATCH -t 01:29:00
#SBATCH -J postprocess
#SBATCH -o slurm.out
#SBATCH -e slurm.err
#SBATCH -A m2705
#SBATCH -N 1
#SBATCH --ntasks-per-node=128
#SBATCH -q regular

module load python
python -u calculate_coh_sqe.py > output.txt
