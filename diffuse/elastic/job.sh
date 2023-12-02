#!/bin/bash -l
#SBATCH -C cpu
#SBATCH -t 00:30:00
#SBATCH -J diffuse
#SBATCH -o slurm.out
#SBATCH -e slurm.err
#SBATCH -A m2705
#SBATCH -N 1
#SBATCH --ntasks-per-node=128
#SBATCH -q debug

#module load python
#conda activate deepmd
#module load openmpi
#srun --mpi=pmi2 /global/homes/x/xinghe/.conda/envs/deepmd/bin/lmp -in in.lammps > log.lammps 

#module load gsl
#module load cray-hdf5-parallel
#module load cray-fftw
#module load lammps
module load python
python -u calculate_elastic_scattering.py > output.txt
