#include <bits/stdc++.h>
#include <mpi.h>
using namespace std;

class Point {
public:
    double x; 
    double y;
};

vector<Point> generate_points(int); 
vector<Point> merge_two_vectors(vector<Point> &A, vector<Point> &B, bool sort_by_x);
vector<Point> parallel_merge_sort(vector<Point> &A, int world_rank, int world_size, bool sort_by_x);
void merge_sort(vector<Point> &A, int start, int end, bool sort_by_x);
void merge(vector<Point> &A, int start, int end, int mid, bool sort_by_x);

int N;

int main(int argc, char **argv) {
    if (argc == 0)
        return 0;
    N = atoi(argv[1]);

    MPI_Init(NULL, NULL);
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
<
    vector<Point> points;
    if (world_rank == 0)
        points = generate_points(N);
  
    vector<Point> points_x = parallel_merge_sort(points, world_rank, world_size, true); // points sorted by x for process 0
    points.clear();

    
    
    MPI_Finalize();
    return 0;
}

vector<Point> parallel_merge_sort(vector<Point> &points, int world_rank, int world_size, bool sort_by_x) {
    // custom data type
    MPI_Datatype dt_point;
    MPI_Type_contiguous(2, MPI_DOUBLE, &dt_point);
    MPI_Type_commit(&dt_point);

    vector<Point> local_points(N / world_size); // local buffer to store a piece of the vector
    MPI_Scatter(points.data(), N / world_size, dt_point, local_points.data(), N / world_size, dt_point, 0, MPI_COMM_WORLD);
    
    merge_sort(local_points, 0, (N / world_size) - 1, sort_by_x);
    
    int size = world_size, rank = world_rank, msg_size = N / world_size;
    while (size != 1) {
        if (rank < size) {
            if (rank >= size / 2) {
                MPI_Send(local_points.data(), msg_size, dt_point, rank - size / 2, 0, MPI_COMM_WORLD);
            } else {
                vector<Point> tmp(msg_size);
                MPI_Recv(tmp.data(), msg_size, dt_point, rank + size / 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                local_points = merge_two_vectors(local_points, tmp, sort_by_x);
            }
        }
        size /= 2;
        msg_size *= 2;
    }

    if (world_rank != 0)
        local_points.clear();

    return local_points;
}

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

void merge_sort(vector<Point> &A, int start, int end, bool sort_by_x) {
    if (start < end) {
        int mid = floor((start + end) / 2);
        merge_sort(A, start, mid, sort_by_x);
        merge_sort(A, mid + 1, end, sort_by_x);
        merge(A, start, end, mid, sort_by_x);
    }
}

void merge(vector<Point> &A, int start, int end, int mid, bool sort_by_x) {
    int i = start, j = mid + 1, k = 0;
    vector<Point> B(end - start + 1);
 
    while (i <= mid && j <= end) {
        if (sort_by_x) {
            if (A[i].x <= A[j].x) 
                B[k++] = A[i++];
            else
                B[k++] = A[j++];
        } else {
            if (A[i].y <= A[j].y) 
                B[k++] = A[i++];
            else
                B[k++] = A[j++];
        }
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