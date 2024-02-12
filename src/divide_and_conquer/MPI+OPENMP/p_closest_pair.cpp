#define _GNU_SOURCE // sched_getcpu(3) is glibc-specific (see the man page)
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sched.h>
#include <mpi.h>
#include <vector>
#include <algorithm>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

class Point {
public:
    double x; 
    double y;
};

int N;

vector<Point> generate_points(int); 
void parallel_merge_sort(vector<Point> &points, int rank, int size, MPI_Datatype dt_point, int threads);
vector<Point> merge_two_vectors(vector<Point> &A, vector<Point> &B, bool sort_by_x);
void merge_two_vectors_in_place(vector<Point> &src, vector<Point> &dst, int i, int j, int k, int n);

int main(int argc, char **argv) {
    if (argc == 1)
        return 0;
    N = atoi(argv[1]);

    // initialize MPI
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // custom data type
    MPI_Datatype dt_point;
    MPI_Type_contiguous(2, MPI_DOUBLE, &dt_point);
    MPI_Type_commit(&dt_point);

    char _hostname[MPI_MAX_PROCESSOR_NAME];
    int _hostname_len;
    MPI_Get_processor_name(_hostname, &_hostname_len);

    int threads = omp_get_max_threads() / size;

    if (rank==0) {
        printf("NUMBER OF PROCESSES %d\n",size);
        printf("THREADS PER PROCESS %d\n",threads);
    }

    // root process holds the data
    vector<Point> points;
    if (rank == 0)
        points = generate_points(N);
    
    parallel_merge_sort(points, rank, size, dt_point, threads);

    if (rank == 0) {
        cout << "POINTS AFTER SORTING" << endl;
        for (int i = 0; i < points.size(); i++) {
            cout << points[i].x << " , " << points[i].y << endl;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

/**
 * @brief Sort a vector of points in parallel by the x-coordinate.
 * 
 * @param points The points to sort
 * @param rank The process rank
 * @param size The number of total processes
 * @param dt_point MPI custom data type to use
*/
void parallel_merge_sort(vector<Point> &points, int rank, int size, MPI_Datatype dt_point, int threads) {
    vector<Point> local_points(N / size);
    MPI_Scatter(points.data(), N / size, dt_point, local_points.data(), N / size, dt_point, 0, MPI_COMM_WORLD);

    #pragma omp parallel num_threads(threads)
    {
        int thread_id = omp_get_thread_num();
        int n_points_to_sort = (N / size) / threads;

        if (thread_id != threads - 1)
            sort(local_points.begin() + thread_id * n_points_to_sort, local_points.begin() + thread_id * n_points_to_sort + n_points_to_sort, 
            [](Point const & a, Point const & b) -> bool
            { return a.x < b.x; } 
            );
        else 
            sort(local_points.begin() + thread_id * n_points_to_sort, local_points.end(), 
            [](Point const & a, Point const & b) -> bool
            { return a.x < b.x; } 
            );
    }
    #pragma omp barrier

    vector<Point> local_points_copy(local_points.size());

    #pragma omp parallel num_threads(threads)
    {
        int thread_id = omp_get_thread_num();
        bool turn = true; // used to decide whether to write in local_points or local_points_copy
        int msg_size = local_points.size() / threads, t_size = threads;

        while (t_size != 1) {
            if (thread_id < t_size && thread_id < t_size / 2) {
                if (!turn) {
                    merge_two_vectors_in_place(local_points_copy, local_points, thread_id * msg_size, (thread_id + t_size / 2) * msg_size, thread_id * msg_size * 2, msg_size);
                } else { 
                    merge_two_vectors_in_place(local_points, local_points_copy, thread_id * msg_size, (thread_id + t_size / 2) * msg_size, thread_id * msg_size * 2, msg_size);
                }
            }
            t_size /= 2;
            msg_size *= 2;
            turn = !turn;
            #pragma omp barrier 
        }
        if (!turn && thread_id == 0)
            local_points = local_points_copy;
    }

    int msg_size = N / size;
    while (size != 1) {
        if (rank < size) {
            if (rank >= size / 2) {
                MPI_Send(local_points.data(), msg_size, dt_point, rank - size / 2, 0, MPI_COMM_WORLD);
            } else {
                vector<Point> tmp(msg_size);
                MPI_Recv(tmp.data(), msg_size, dt_point, rank + size / 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                local_points = merge_two_vectors(local_points, tmp, true);
            }
        }
        size /= 2;
        msg_size *= 2;
    }

    if (rank == 0)
        points = local_points;
}

/**
 * @brief Generate n randon points
 * 
 * @param n The number of points to generate
 * @return a vector containing n random points
*/
vector<Point> generate_points(int n) {
    vector<Point> points(n);

    srand(time(NULL));
    const long max_rand = 100000000L;
    double lower_bound = -1000000.0;
    double upper_bound = 1000000.0;

    for (int i = 0; i < n; i++) {
        Point p;
        p.x = lower_bound + (upper_bound - lower_bound) * (rand() % max_rand) / max_rand;
        p.y = lower_bound + (upper_bound - lower_bound) * (rand() % max_rand) / max_rand;
        points[i] = p;
    }
    return points;
}
/**
 * @brief Merge two vectors of size n starting at index i and j from src to dst from index k to k + 2n
 * 
 * @param src The vector to copy from
 * @param dst The vector to copy to
 * @param i Index of the first vector in src
 * @param j Index of the second vector in src
 * @param k Index where to start writing in dst
 * @param n Lenght of the two vectors to merge
*/
void merge_two_vectors_in_place(vector<Point> &src, vector<Point> &dst, int i, int j, int k, int n) {
    int i_max = i + n, j_max = j + n;

    while (i < i_max && j < j_max) {
        if (src[i].x < src[j].x)
            dst[k++] = src[i++];
        else
            dst[k++] = src[j++];
    }

    if (i == i_max) {
        while (j < j_max)
            dst[k++] = src[j++];
    } else {
        while(i < i_max)    
            dst[k++] = src[i++];
    }
}

/**
 * @brief Merge vectors A and B
 * 
 * @param A The first vector
 * @param B The second vector
 * @param sort_by_x Indicates whether to sort by x or y
 * 
 * @return Vectors A and B merged together.
*/
vector<Point> merge_two_vectors(vector<Point> &A, vector<Point> &B, bool sort_by_x) {
    int n = A.size();
    vector<Point> C(n * 2);
    int i = 0, j = 0, k = 0;

    while (i < n && j < n) {
        if (sort_by_x) {
            if (A[i].x <= B[j].x)
                C[k++] = A[i++];
            else    
                C[k++] = B[j++];
        } else {
            if (A[i].y <= B[j].y)
                C[k++] = A[i++];
            else    
                C[k++] = B[j++];
        }
    }

    if (i == n)
        while (j < n)
            C[k++] = B[j++];
    else 
        while (i < n)
            C[k++] = A[i++];

    return C;
}