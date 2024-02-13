#!/bin/bash

#PBS -l select=1:ncpus=4:mem=4gb

#set max execution time
#PBS -l walltime=0:15:00

#PBS -q short_cpuQ

mpiexec -n 4 ./openmp_and_mpi_closest_pair 2097152