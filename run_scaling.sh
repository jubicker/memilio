#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 56
#SBATCH -t 6-0:00:00
#SBATCH --exclusive
#SBATCH --output=%A_scaling.out
#SBATCH --error=%A_scaling.err
#SBATCH --nodelist="be-cpu02, be-cpu03, be-cpu04"

for i in 56 32 16 8 4 2 1
do
    echo Run with $i cores.
    srun --cpu-bind=core --cpus-per-task=$i ./build/bin/$1
done
