#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 128
#SBATCH -t 0-24:00:00
#SBATCH --exclusive
#SBATCH --output=%A.out
#SBATCH --error=%A.err
#SBATCH --account=loki
#SBATCH --threads-per-core=2

for i in 1 2 4 8 16 32 64 128
do
    echo Run with $i cores.
    for j in {0..0}
    do
        srun --cpu-bind=core --cpus-per-task=$i ./cpp/build/bin/$1
    done
done
