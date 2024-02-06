#!/bin/bash

#PBS -l select=1:ncpus=1:mem=4gb

#set max execution time
#PBS -l walltime=0:15:00

#imposta la coda di esecuzione
#PBS -q short_cpuQ

module load mpich-3.2
#mpic++ -g -Wall -o ./broadcast/broadcast ./broadcast/broadcast.c
#mpirun.actual -n 2 ./broadcast/broadcast
./project/naive