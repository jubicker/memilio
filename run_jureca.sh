#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 50
#SBATCH -t 0-24:00:00
#SBATCH --exclusive
#SBATCH --output=%A.out
#SBATCH --error=%A.err
#SBATCH --account=loki
#SBATCH --threads-per-core=1
srun --cpu-bind=core --cpus-per-task=50 ./build/bin/$1
