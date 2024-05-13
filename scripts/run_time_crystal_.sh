#!/bin/bash

mkdir nohup

TASK_SET=002

n_tasks=1
n_cpus_task=4

dt=0.01
T=15
avg=4000

module purge > /dev/null 2>&1

eval "$(conda shell.bash hook)"
conda activate dev_phys

for TASK_SET in 004 005
do
	#MIXED PHASE
	TASKNAME="MIX"

	linear=0
	M=1
	v0=1
	eta_1=1
	eta_2=1
	T1=1
	T2=1

	for alpha in 0.05 0.1 0.3 0.5 0.7 0.9
	do
	  out_name=set$TASK_SET-avg$avg-dt$dt-T$T-linear$linear-eta1$eta_1-eta2$eta_2-T1$T1-T2$T2-v0$v0-M$M-alpha$alpha
	  JOBNAME=$TASKNAME$alpha
	  srun -n $n_tasks -c $n_cpus_task --time=100:00:00 --job-name=$JOBNAME python scripts/msd_fle_time_crystal_parallel.py $alpha $linear $M $v0 $eta_1 $eta_2 $T1 $T2 $dt $T $avg $TASK_SET > nohup/$out_name.out 2>&1 &
	done

	#IN EQUILIBRIUM
	TASKNAME="EQ"

	linear=1
	M=1
	v0=1
	eta_1=1
	eta_2=1
	T1=1
	T2=1

	for alpha in 0.05 0.1 0.3 0.5 0.7 0.9
	do
	  out_name=set$TASK_SET-avg$avg-dt$dt-T$T-linear$linear-eta1$eta_1-eta2$eta_2-T1$T1-T2$T2-v0$v0-M$M-alpha$alpha
	  JOBNAME=$TASKNAME$alpha
	  srun -n $n_tasks -c $n_cpus_task --time=100:00:00 --job-name=$JOBNAME python scripts/msd_fle_time_crystal_parallel.py $alpha $linear $M $v0 $eta_1 $eta_2 $T1 $T2 $dt $T $avg $TASK_SET > nohup/$out_name.out 2>&1 &
	done
done
 
wait
