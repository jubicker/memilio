#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 128
#SBATCH -t 0-24:00:00
#SBATCH --exclusive
#SBATCH --output=%A_study.out
#SBATCH --error=%A_study.err
#SBATCH --account=loki
#SBATCH --threads-per-core=2

srun --cpu-bind=core --cpus-per-task=128 ./build/bin/$1
