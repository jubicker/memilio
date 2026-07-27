#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 128
#SBATCH -t 1-00:00:00
#SBATCH --exclusive
#SBATCH --output=%A_fitting_ag.out
#SBATCH --error=%A_fitting_ag.err
#SBATCH --account=loki
#SBATCH --threads-per-core=1

module load Python/3.13.5

cd /p/project1/loki/bicker1/memilio/virt_env_ag/bin
. ./activate
cd /p/project1/loki/bicker1/memilio

# pyabc's default sampler picks its worker-process count from os.cpu_count(),
# which reports this node's full 256 logical threads and ignores the SLURM
# allocation. Pin it to the cores actually granted to this job.
export PYABC_NUM_PROCS=$SLURM_CPUS_PER_TASK
# The simulation binary is built with OpenMP. Without this, each of the
# PYABC_NUM_PROCS worker processes would also spawn its own OpenMP thread
# team, oversubscribing the allocated cores many times over.
export OMP_NUM_THREADS=1

srun --cpu-bind=core --cpus-per-task=128 python /p/project1/loki/bicker1/memilio/cpp/simulations/hybrid_simulations/sir_metapop/bindings/run_fitting_germany_age_groups.py
