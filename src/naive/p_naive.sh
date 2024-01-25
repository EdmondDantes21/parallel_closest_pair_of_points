#!/bin/bash

#PBS -l select=1:ncpus=64:mem=4gb

#set max execution time
#PBS -l walltime=0:18:00

#imposta la coda di esecuzione
#PBS -q short_cpuQ

module load mpich-3.2

printf "P = 1, N = 1024" && mpirun.actual -n 1 ./project/p_naive 1024 &&
printf "P = 1, N = 2048" && mpirun.actual -n 1 ./project/p_naive 2048 &&
printf "P = 1, N = 4096" && mpirun.actual -n 1 ./project/p_naive 4096 &&
printf "P = 1, N = 8192" && mpirun.actual -n 1 ./project/p_naive 8192 &&
printf "P = 1, N = 16384" && mpirun.actual -n 1 ./project/p_naive 16384 &&
printf "P = 1, N = 32768" && mpirun.actual -n 1 ./project/p_naive 32768 &&
printf "P = 1, N = 65536" && mpirun.actual -n 1 ./project/p_naive 65536 &&
printf "P = 1, N = 131072" && mpirun.actual -n 1 ./project/p_naive 131072 &&
printf "P = 1, N = 262144" && mpirun.actual -n 1 ./project/p_naive 262144 &&
printf "P = 2, N = 1024" && mpirun.actual -n 2 ./project/p_naive 1024 &&
printf "P = 2, N = 2048" && mpirun.actual -n 2 ./project/p_naive 2048 &&
printf "P = 2, N = 4096" && mpirun.actual -n 2 ./project/p_naive 4096 &&
printf "P = 2, N = 8192" && mpirun.actual -n 2 ./project/p_naive 8192 &&
printf "P = 2, N = 16384" && mpirun.actual -n 2 ./project/p_naive 16384 &&
printf "P = 2, N = 32768" && mpirun.actual -n 2 ./project/p_naive 32768 &&
printf "P = 2, N = 65536" && mpirun.actual -n 2 ./project/p_naive 65536 &&
printf "P = 2, N = 131072" && mpirun.actual -n 2 ./project/p_naive 131072 &&
printf "P = 2, N = 262144" && mpirun.actual -n 2 ./project/p_naive 262144 &&
printf "P = 4, N = 1024" && mpirun.actual -n 4 ./project/p_naive 1024 &&
printf "P = 4, N = 2048" && mpirun.actual -n 4 ./project/p_naive 2048 &&
printf "P = 4, N = 4096" && mpirun.actual -n 4 ./project/p_naive 4096 &&
printf "P = 4, N = 8192" && mpirun.actual -n 4 ./project/p_naive 8192 &&
printf "P = 4, N = 16384" && mpirun.actual -n 4 ./project/p_naive 16384 &&
printf "P = 4, N = 32768" && mpirun.actual -n 4 ./project/p_naive 32768 &&
printf "P = 4, N = 65536" && mpirun.actual -n 4 ./project/p_naive 65536 &&
printf "P = 4, N = 131072" && mpirun.actual -n 4 ./project/p_naive 131072 &&
printf "P = 4, N = 262144" && mpirun.actual -n 4 ./project/p_naive 262144 &&
printf "P = 8, N = 1024" && mpirun.actual -n 8 ./project/p_naive 1024 &&
printf "P = 8, N = 2048" && mpirun.actual -n 8 ./project/p_naive 2048 &&
printf "P = 8, N = 4096" && mpirun.actual -n 8 ./project/p_naive 4096 &&
printf "P = 8, N = 8192" && mpirun.actual -n 8 ./project/p_naive 8192 &&
printf "P = 8, N = 16384" && mpirun.actual -n 8 ./project/p_naive 16384 &&
printf "P = 8, N = 32768" && mpirun.actual -n 8 ./project/p_naive 32768 &&
printf "P = 8, N = 65536" && mpirun.actual -n 8 ./project/p_naive 65536 &&
printf "P = 8, N = 131072" && mpirun.actual -n 8 ./project/p_naive 131072 &&
printf "P = 8, N = 262144" && mpirun.actual -n 8 ./project/p_naive 262144 &&
printf "P = 16, N = 1024" && mpirun.actual -n 16 ./project/p_naive 1024 &&
printf "P = 16, N = 2048" && mpirun.actual -n 16 ./project/p_naive 2048 &&
printf "P = 16, N = 4096" && mpirun.actual -n 16 ./project/p_naive 4096 &&
printf "P = 16, N = 8192" && mpirun.actual -n 16 ./project/p_naive 8192 &&
printf "P = 16, N = 16384" && mpirun.actual -n 16 ./project/p_naive 16384 &&
printf "P = 16, N = 32768" && mpirun.actual -n 16 ./project/p_naive 32768 &&
printf "P = 16, N = 65536" && mpirun.actual -n 16 ./project/p_naive 65536 &&
printf "P = 16, N = 131072" && mpirun.actual -n 16 ./project/p_naive 131072 &&
printf "P = 16, N = 262144" && mpirun.actual -n 16 ./project/p_naive 262144 &&
printf "P = 32, N = 1024" && mpirun.actual -n 32 ./project/p_naive 1024 &&
printf "P = 32, N = 2048" && mpirun.actual -n 32 ./project/p_naive 2048 &&
printf "P = 32, N = 4096" && mpirun.actual -n 32 ./project/p_naive 4096 &&
printf "P = 32, N = 8192" && mpirun.actual -n 32 ./project/p_naive 8192 &&
printf "P = 32, N = 16384" && mpirun.actual -n 32 ./project/p_naive 16384 &&
printf "P = 32, N = 32768" && mpirun.actual -n 32 ./project/p_naive 32768 &&
printf "P = 32, N = 65536" && mpirun.actual -n 32 ./project/p_naive 65536 &&
printf "P = 32, N = 131072" && mpirun.actual -n 32 ./project/p_naive 131072 &&
printf "P = 32, N = 262144" && mpirun.actual -n 32 ./project/p_naive 262144 &&
printf "P = 64, N = 1024" && mpirun.actual -n 64 ./project/p_naive 1024 &&
printf "P = 64, N = 2048" && mpirun.actual -n 64 ./project/p_naive 2048 &&
printf "P = 64, N = 4096" && mpirun.actual -n 64 ./project/p_naive 4096 &&
printf "P = 64, N = 8192" && mpirun.actual -n 64 ./project/p_naive 8192 &&
printf "P = 64, N = 16384" && mpirun.actual -n 64 ./project/p_naive 16384 &&
printf "P = 64, N = 32768" && mpirun.actual -n 64 ./project/p_naive 32768 &&
printf "P = 64, N = 65536" && mpirun.actual -n 64 ./project/p_naive 65536 &&
printf "P = 64, N = 131072" && mpirun.actual -n 64 ./project/p_naive 131072 &&
printf "P = 64, N = 262144" && mpirun.actual -n 64 ./project/p_naive 262144


