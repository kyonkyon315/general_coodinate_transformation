#!/bin/bash
#SBATCH -p gr20001a
#SBATCH -t 00:30:00
#SBATCH -o py_out.txt
#SBATCH -e py_err.txt
#SBATCH --rsc p=1:t=1:c=32:m=4G

set -e

. /usr/share/Modules/init/bash

# module purge しない（重要）

# 自分のvenvを有効化
source ~/pyenv/bin/activate

echo "Python:"
which python
python --version

echo "SLURM_CPUS_PER_TASK = $SLURM_CPUS_PER_TASK"
echo "SLURM_NTASKS = $SLURM_NTASKS"
echo "SLURM_JOB_CPUS_PER_NODE = $SLURM_JOB_CPUS_PER_NODE"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

srun -n 1 python two_stream_instability_super.py
