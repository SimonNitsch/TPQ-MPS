#!/bin/sh
#PBS -l nodes=1:ppn=16
#PBS -j oe
cd $PBS_O_WORKDIR
NPROCS=1

export I_MPI_FABRICS=tcp
MPIEXEC=mpirun
MPI_PROGRAM=tdvp

mpirun -machinefile $PBS_NODEFILE -np $NPROCS ./tdvp
