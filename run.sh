#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 56
#SBATCH -t 6-0:00:00
#SBATCH --exclusive
#SBATCH --output=smm-%A.out
#SBATCH --error=smm-%A.err
#SBATCH --nodelist="be-cpu02, be-cpu03, be-cpu04"
srun --cpu-bind=core --cpus-per-task=56 ./build/bin/$1
