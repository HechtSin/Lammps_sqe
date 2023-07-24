#!/bin/bash -l
#SBATCH -C cpu
#SBATCH -t 00:29:00
#SBATCH -J LAMMPS
#SBATCH -o slurm.out
#SBATCH -e slurm.err
#SBATCH -A m2705
#SBATCH -N 1
#SBATCH --ntasks-per-node=128
#SBATCH -q debug

module load gsl
module load cray-hdf5-parallel
module load cray-fftw
module load lammps
srun -n 128 -c 2 --cpu_bind=cores lmp -in in.lammps > log.lammps 
