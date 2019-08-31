#!/bin/bash
#SBATCH --time=0-240:0
#SBATCH --cpus-per-task=16
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
./fluxos_cpp