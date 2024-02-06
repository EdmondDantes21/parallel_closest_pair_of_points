#include <bits/stdc++.h>
#include <iostream>
#include <mpi.h>
using namespace std;

class Point {
public:
    double x; 
    double y;

    bool operator <=(const Point& b) {
        return x <= b.x;
    }
};

vector<Point> generate_points(int);
void merge_sort(vector<Point> &A, int start, int end);
void merge(vector<Point> &A, int start, int end, int mid);
vector<Point> merge_two_vectors(vector<Point> &A, vector<Point> &B);    

int main(int argc, char** argv) {
    if (argc == 0)
        return 0;
    int N = atoi(argv[1]);

    MPI_Init(NULL, NULL);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // custom data type
    MPI_Datatype dt_point;
    MPI_Type_contiguous(2, MPI_DOUBLE, &dt_point);
    MPI_Type_commit(&dt_point);

    vector<Point> points; 
    if (world_rank == 0)    // root process generates data
        points = generate_points(N);

    vector<Point> local_points(N / world_size); // local buffer to store a piece of the vector
    MPI_Scatter(points.data(), N / world_size, dt_point, local_points.data(), N / world_size, dt_point, 0, MPI_COMM_WORLD);
    
    merge_sort(local_points, 0, (N / world_size) - 1);
    
    int size = world_size, rank = world_rank, msg_size = N / world_size;
    while (size != 1) {
        if (rank < size) {
            if (rank >= size / 2) {
                MPI_Send(local_points.data(), msg_size, dt_point, rank - size / 2, 0, MPI_COMM_WORLD);
            } else {
                vector<Point> tmp(msg_size);
                MPI_Recv(tmp.data(), msg_size, dt_point, rank + size / 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                local_points = merge_two_vectors(local_points, tmp);
            }
        }
        size /= 2;
        msg_size *= 2;
    }
    
    MPI_Finalize();
    return 0;
}

vector<Point> merge_two_vectors(vector<Point> &A, vector<Point> &B) {
    int n = A.size();
    vector<Point> C(n * 2);
    int i = 0, j = 0, k = 0;

    while (i < n && j < n) {
        if (A[i] <= B[j])
            C[k++] = A[i++];
        else    
            C[k++] = B[j++];
    }

    if (i == n)
        while (j < n)
            C[k++] = B[j++];
    else 
        while (i < n)
            C[k++] = A[i++];

    return C;
}

void merge_sort(vector<Point> &A, int start, int end) {
    if (start < end) {
        int mid = floor((start + end) / 2);
        merge_sort(A, start, mid);
        merge_sort(A, mid + 1, end);
        merge(A, start, end, mid);
    }
}

void merge(vector<Point> &A, int start, int end, int mid) {
    int i = start, j = mid + 1, k = 0;
    vector<Point> B(end - start + 1);
 
    while (i <= mid && j <= end) {
        if (A[i] <= A[j]) 
            B[k++] = A[i++];
        else
            B[k++] = A[j++];
    }

    if (i > mid) 
        while (j <= end)
            B[k++] = A[j++];
    else 
        while (i <= mid)
            B[k++] = A[i++];

    for (int l = 0; l < k; l++)
        A[start++] = B[l];
}

vector<Point> generate_points(int n) {
    vector<Point> points(n);

    srand(time(NULL));
    const long max_rand = 1000000L;
    double lower_bound = -10000.0;
    double upper_bound = 10000.0;

    for (int i = 0; i < n; i++) {
        Point p;
        p.x = lower_bound + (upper_bound - lower_bound) * (rand() % max_rand) / max_rand;
        p.y = lower_bound + (upper_bound - lower_bound) * (rand() % max_rand) / max_rand;
        points[i] = p;
    }
    return points;
}