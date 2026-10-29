#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 128
#SBATCH -t 0-1:00:00
#SBATCH --exclusive
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --account=loki

echo $1

for i in 32 #1 2 4 8 16 32 64 128
do
    for j in {0..4}
    do
        echo Run with $i cores.
        srun --cpu-bind=core --cpus-per-task=$i ./build/bin/$1
    done
done
