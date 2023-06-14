#!/bin/bash
#------- qsub option -----------
#PBS -q DBG
#PBS -N calc
#PBS --group=K2319
#PBS -l elapstim_req=00:10:00
#PBS -l cpunum_job=4
#PBS -v OMP_NUM_THREADS=4
#PBS -T intmpi
#PBS -b 1

export LANG="en_US.utf8"
export PYTHONUNBUFFERED=1
cd $PBS_O_WORKDIR
exec 1<&- && exec 2<&-
exec 1<>${PBS_JOBNAME}_LOGFILE.log && exec 2>&1

# ------- Program execution -----------
module load BaseCPU/2023
source activate catlearn

echo "--------------ENVIRONMENT------------------"
echo "LANGUAGE:  " $LANG
echo "OMP_NUM_THREADS:  " $OMP_NUM_THREADS
echo "CONDA ENV:  " $CONDA_DEFAULT_ENV
echo "PYTHON:  " $(which python)
echo "PYTHON VERSION:  " $(python -V -V)
echo "LOCATION:  " $(pwd)
echo "-------------------------------------------"

echo "========= Job started  at `date` =========="
SCRIPT_FILE="catlearn-example-01.py"
mkdir -p run_dir && cd run_dir
cp ../*traj .
time python ../$SCRIPT_FILE

echo "========= Job finished at `date` =========="