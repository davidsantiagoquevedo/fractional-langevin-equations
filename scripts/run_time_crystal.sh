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

#TIME GLASS
TASKNAME="TG"

linear="False"
M=1
v0=1
eta_1=1
eta_2=0
T1=1
T2=0

for alpha in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.3 0.5 0.7 0.9 0.99
do
  out_name=set$TASK_SET-avg$avg-dt$dt-T$T-linear$linear-eta1$eta_1-eta2$eta_2-T1$T1-T2$T2-v0$v0-M$M-alpha$alpha
  JOBNAME=$TASKNAME$alpha
  srun -n $n_tasks -c $n_cpus_task --time=100:00:00 --job-name=$JOBNAME python scripts/msd_fle_time_crystal_parallel.py $alpha $linear $M $v0 $dt $T $avg $TASK_SET > nohup/$out_name.out 2>&1 &
done

#TIME CRYSTAL
TASKNAME="TC"

linear="False"
M=1
v0=1
eta_1=0
eta_2=1
T1=0
T2=1

for alpha in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.3 0.5 0.7 0.9 0.99
do
  out_name=set$TASK_SET-avg$avg-dt$dt-T$T-linear$linear-eta1$eta_1-eta2$eta_2-T1$T1-T2$T2-v0$v0-M$M-alpha$alpha
  JOBNAME=$TASKNAME$alpha
  srun -n $n_tasks -c $n_cpus_task --time=100:00:00 --job-name=$JOBNAME python scripts/msd_fle_time_crystal_parallel.py $alpha $linear $M $v0 $dt $T $avg $TASK_SET > nohup/$out_name.out 2>&1 &
done

#TIME CRYSTAL + GLASSES
TASKNAME="TCG"

linear="False"
M=1
v0=1
eta_1=1
eta_2=1
T1=1
T2=1

for alpha in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.3 0.5 0.7 0.9 0.99
do
  out_name=set$TASK_SET-avg$avg-dt$dt-T$T-linear$linear-eta1$eta_1-eta2$eta_2-T1$T1-T2$T2-v0$v0-M$M-alpha$alpha
  JOBNAME=$TASKNAME$alpha
  srun -n $n_tasks -c $n_cpus_task --time=100:00:00 --job-name=$JOBNAME python scripts/msd_fle_time_crystal_parallel.py $alpha $linear $M $v0 $dt $T $avg $TASK_SET > nohup/$out_name.out 2>&1 &
done

#IN EQUILIBRIUM
TASKNAME="EQ"

linear="True"
M=1
v0=1
eta_1=1
eta_2=1
T1=1
T2=1

for alpha in 0.05 0.1 0.3 0.5 0.7 0.9 0.99
do
  out_name=set$TASK_SET-avg$avg-dt$dt-T$T-linear$linear-eta1$eta_1-eta2$eta_2-T1$T1-T2$T2-v0$v0-M$M-alpha$alpha
  JOBNAME=$TASKNAME$alpha
  srun -n $n_tasks -c $n_cpus_task --time=100:00:00 --job-name=$JOBNAME python scripts/msd_fle_time_crystal_parallel.py $alpha $linear $M $v0 $dt $T $avg $TASK_SET > nohup/$out_name.out 2>&1 &
done
 
wait