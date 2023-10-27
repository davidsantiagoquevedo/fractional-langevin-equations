#!/bin/bash

module purge > /dev/null 2>&1

eval "$(conda shell.bash hook)"
conda activate dev_phys

srun srun -n 1 -c 4 python scripts/msd_fle_two_baths_parallel.py 0.95 1 1 0 1 0 > nohup0.out 2>&1 &
srun srun -n 1 -c 4 python scripts/msd_fle_two_baths_parallel.py 0.955 1 1 0 1 0 > nohup0.out 2>&1 &

wait