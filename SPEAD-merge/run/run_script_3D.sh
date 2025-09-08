#!/bin/sh
#SBATCH --job-name=SPEAD_global_cont
#SBATCH --time=04-12:00:00     # run no longer than 60 hours
#SBATCH --mem-per-cpu=2500
#SBATCH --ntasks=24
#SBATCH --output=/home/legland/LEGLAND/SPEAD-merge/output/meta-data/SLURM_OUTPUT
#SBATCH --error=/home/legland/LEGLAND/SPEAD-merge/output/meta-data/SLURM_ERROR
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillaume.le-gland@univ-amu.fr


export initdir=/home/legland/LEGLAND/SPEAD-merge/run/
export rundir=/home/legland/LEGLAND/SPEAD-merge/output/meta-data/darwin_latest/
mkdir $rundir
cd $rundir
rm *

cp $initdir/mitgcmuv $rundir
cp $initdir/POLY3.COEFFS $rundir
cp $initdir/../data/* $rundir
cp $initdir/pickup* $rundir
sleep 10s

mpirun -np 24 ./mitgcmuv
sleep 10s

