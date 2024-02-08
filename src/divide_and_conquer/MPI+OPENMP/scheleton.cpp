#define _GNU_SOURCE // sched_getcpu(3) is glibc-specific (see the man page)
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sched.h>
#include <mpi.h>

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char **argv){
    // initialize MPI
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    char _hostname[MPI_MAX_PROCESSOR_NAME];
    int _hostname_len;
    MPI_Get_processor_name(_hostname, &_hostname_len);

    int threads = omp_get_max_threads();

    if (rank==0) {
        printf("Processes for MPI %d\n",size);
        printf("Threading for OMP %d\n",threads);
    }

    #pragma omp barrier
    #pragma omp parallel num_threads(threads)
    {
        int thread_id = omp_get_thread_num();
        int cpu_num = sched_getcpu();
        printf("Node %s, Thread %d of process %d run on CPU %d \n", _hostname, thread_id, rank, cpu_num);
    }
    #pragma omp barrier

    sleep(1);
    printf("Process rank %d is done \n", rank);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}