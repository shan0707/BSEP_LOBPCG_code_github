#!/bin/bash
#SBATCH --job-name=test_sparse         # 作业名称
#SBATCH --cpus-per-task=64 # 每个任务申请的 CPU 核心数
#SBATCH --time=12:00:00           # 作业运行时间
#SBATCH --nodelist=bigMem1

module load MATLAB
matlab -nodesktop -nosplash -r test_sparse