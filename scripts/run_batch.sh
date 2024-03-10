#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=10:00:00 

module purge > /dev/null 2>&1

eval "$(conda shell.bash hook)"
conda activate dev_phys


srun --exclusive -n 1 python scripts/msd_fle_two_baths_parallel.py 0.95 1 1 0 1 0 > nohup1.out 2>&1 & 
srun --exclusive -n 1 python scripts/msd_fle_two_baths_parallel.py 0.955 1 1 0 1 0 > nohup2.out 2>&1 & 
srun --exclusive -n 1 python scripts/msd_fle_two_baths_parallel.py 0.96 1 1 0 1 0 > nohup3.out 2>&1 & 
srun --exclusive -n 1 python scripts/msd_fle_two_baths_parallel.py 0.965 1 1 0 1 0 > nohup4.out 2>&1 & 
#srun --exclusive -n 1 python scripts/msd_fle_two_baths_parallel.py 0.97 1 1 0 1 0 > nohup5.out 2>&1 & 
#srun --exclusive -n 1 python scripts/msd_fle_two_baths_parallel.py 0.975 1 1 0 1 0 > nohup6.out 2>&1 & 

wait

