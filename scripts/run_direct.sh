#!/bin/bash

mkdir nohup

TASK_SET=001

n_tasks=1
n_cpus_task=4

dt=0.05
T=100
avg=4000

module purge > /dev/null 2>&1

eval "$(conda shell.bash hook)"
conda activate dev_phys

# TIME GLASS
# TASKNAME="TG"
# 
# A=1.0
# eta=10.0
# C=0.0
# theta_H=0.0
# theta_12=1.0

# for H in 0.95 0.952 0.954 0.956 0.958 0.96 0.962 0.964 0.966 0.968 0.97 0.972 0.974 
# do
#   out_name=$TASKNAME-avg$avg-dt$dt-T$T-H$H-A$A-eta$eta-C$C-theta12_$theta_12-thetaH$theta_H
#   JOBNAME=$TASKNAME$H
#   srun -n $n_tasks -c $n_cpus_task --time=100:00:00 --job-name=$JOBNAME python scripts/msd_fle_two_baths_parallel.py $H $A $eta $C $theta_12 $theta_H $dt $T $avg $TASK_SET > nohup/$out_name.out 2>&1 &
# done

# LUTZ
TASKNAME="LU"
# 
A=1.0
eta=1.0
C=0.0
theta_H=1.0
theta_12=0.0

for H in 0.500000001
do
  out_name=$TASKNAME$TASK_SET-avg$avg-dt$dt-T$T-H$H-A$A-eta$eta-C$C-theta12_$theta_12-thetaH$theta_H
  JOBNAME=$TASKNAME$H
  srun -n $n_tasks -c $n_cpus_task --time=100:00:00 --job-name=$JOBNAME python scripts/msd_fle_two_baths_parallel.py $H $A $eta $C $theta_12 $theta_H $dt $T $avg $TASK_SET > nohup/$out_name.out 2>&1 &
done

# oBm + colored
# TASKNAME="OC"

# A=1.0
# eta=0.0
# C=1.0
# theta_H=1.0
# theta_12=0.0

# for H in 0.50001 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 0.99 1
# do
#   out_name=$TASKNAME-avg$avg-dt$dt-T$T-H$H-A$A-eta$eta-C$C-theta12_$theta_12-thetaH$theta_H
#   JOBNAME=$TASKNAME$H
#   srun -n $n_tasks -c $n_cpus_task --time=100:00:00 --job-name=$JOBNAME python scripts/msd_fle_two_baths_parallel.py $H $A $eta $C $theta_12 $theta_H $dt $T $avg $TASK_SET > nohup/$out_name.out 2>&1 &
# done

# White + Coloured
# TASKNAME="WC"

# A=1.0
# eta=1.0
# C=0.0
# theta_H=1.0
# theta_12=1.0

# for H in 0.50001 0.6 0.64 0.68 0.74 0.76 0.8 0.86 0.88 0.9 0.92 0.94
# do
#   out_name=$TASKNAME-avg$avg-dt$dt-T$T-H$H-A$A-eta$eta-C$C-theta12_$theta_12-thetaH$theta_H
#   JOBNAME=$TASKNAME$H
#   srun -n $n_tasks -c $n_cpus_task --time=100:00:00 --job-name=$JOBNAME python scripts/msd_fle_two_baths_parallel.py $H $A $eta $C $theta_12 $theta_H $dt $T $avg $TASK_SET > nohup/$out_name.out 2>&1 &
# done

# All terms
# TASKNAME="AL"
#  
# A=1.0
# eta=1.0
# C=1.0
# theta_H=1.0
# theta_12=1.0
# 
# for H in 0.50001 0.52 0.54 0.56 0.58 0.6 0.62 0.64 0.66 0.68 0.7 0.72 0.74 0.76 0.78 0.8 0.82 0.84 0.86 0.88 0.9 0.92 0.94 0.96 0.98 1
# do
  # out_name=$TASKNAME-avg$avg-dt$dt-T$T-H$H-A$A-eta$eta-C$C-theta12_$theta_12-thetaH$theta_H
  # JOBNAME=$TASKNAME$H
  # srun -n $n_tasks -c $n_cpus_task --time=100:00:00 --job-name=$JOBNAME python scripts/msd_fle_two_baths_parallel.py $H $A $eta $C $theta_12 $theta_H $dt $T $avg $TASK_SET > nohup/$out_name.out 2>&1 &
# done
#  
# wait

# External Damping
# TASKNAME="ED"
 
# A=1.0
# eta=1.0
# C=1.0
# theta_H=1.0
# theta_12=0.0

# for H in 0.50001 0.52 0.54 0.56 0.58 0.6 0.62 0.64 0.66 0.68 0.7 0.72 0.74 0.76 0.78 0.8 0.82 0.84 0.86 0.88 0.9 0.92 0.94 0.96 0.98 1
# do
#   out_name=$TASKNAME-avg$avg-dt$dt-T$T-H$H-A$A-eta$eta-C$C-theta12_$theta_12-thetaH$theta_H
#   JOBNAME=$TASKNAME$H
#   srun -n $n_tasks -c $n_cpus_task --time=100:00:00 --job-name=$JOBNAME python scripts/msd_fle_two_baths_parallel.py $H $A $eta $C $theta_12 $theta_H $dt $T $avg $TASK_SET > nohup/$out_name.out 2>&1 &
# done
 
wait