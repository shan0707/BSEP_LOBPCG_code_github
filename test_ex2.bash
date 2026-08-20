#!/bin/bash
#SBATCH --job-name=test_ex2        # 作业名称
#SBATCH --cpus-per-task=54     # 每个任务申请的 CPU 核心数
#SBATCH --time=5:00:00             # 作业运行时间
#SBATCH --nodelist=bigMem0

module load MATLAB
matlab -nodesktop -nosplash -r test_ex2