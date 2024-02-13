#!/bin/bash

#PBS -l select=2:ncpus=16:mem=4gb

#set max execution time
#PBS -l walltime=0:15:00

#PBS -q short_cpuQ

mpiexec -n 32 ./pure_mpi_closest_pair 8388608