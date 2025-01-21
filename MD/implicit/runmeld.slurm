#!/bin/bash
#SBATCH --job-name=new_md
#SBATCH --partition gpu_long
#SBATCH --nodes=2
#SBATCH --gres=gpu:2080:10
#SBATCH --ntasks=20
#SBATCH --output=job_logs/slurm.out
#SBATCH --error=job_logs/slurm.err
#SBATCH --time=48:00:00

source ~/.bashrc
conda activate meld

module load cuda12.2/toolkit/12.2.2

#module load openmpi/gcc/64/4.1.5
#module load gcc/8.3.0
export PATH=/cm/shared/apps/cuda12.2/toolkit/12.2.2/bin:$PATH
#export LD_LIBRARY_PATH=/data01/software/cuda11.2/lib64/:$LD_LIBRARY_PATH
export OPENMM_CUDA_COMPILER=$(which nvcc)

if [ -f remd.log ] ; then
    echo "remd.log found, restarting the simulation"
    prepare_restart --prepare-run
fi

which launch_remd

mpirun -np 20 launch_remd --debug --platform CUDA
