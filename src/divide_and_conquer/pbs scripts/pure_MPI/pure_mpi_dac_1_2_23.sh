#!/bin/bash

#PBS -l select=1:ncpus=2:mem=4gb

#set max execution time
#PBS -l walltime=0:15:00

#PBS -q short_cpuQ

mpiexec -n 2 ./pure_mpi_closest_pair 8388608