#define _GNU_SOURCE // sched_getcpu(3) is glibc-specific (see the man page)
#include <bits/stdc++.h>
#include <sched.h>
#include <mpi.h>

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
void merge_two_vectors_in_place(vector<Point> &src, vector<Point> &dst, int i, int j, int k, int n, bool sort_by_x);
void sort_using_threads(vector<Point> &A, int threads);
vector<Point> slice_vec(vector<Point> &A, int start, int end);
double min_dist(vector<Point>& x_sorted, vector<Point>& y_sorted);
vector<Point> compute_strip(vector<Point>& y_sorted, double delta, double middle);
double brute_force(vector<Point>& points);
inline double dist(Point& p1, Point& p2);

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

    int threads = omp_get_max_threads() / size;

    // root process holds the data
    vector<Point> points;
    if (rank == 0)
        points = generate_points(N);
    
    parallel_merge_sort(points, rank, size, dt_point, threads);

    vector<Point> local_x(N / size);
    MPI_Scatter(points.data(), N / size, dt_point, local_x.data(), N / size, dt_point, 0, MPI_COMM_WORLD);
    vector<Point> local_y = local_x;
    sort_using_threads(local_y, threads);

    if (rank == 0)
        points.clear();

    double local_delta = DBL_MAX;
    {
        vector<double> delta(threads);
        vector<vector<Point>> strips(threads);
        #pragma omp parallel num_threads(threads) reduction(min:local_delta)
        {
            int thread_id = omp_get_thread_num();
            int points_per_thread = local_x.size() / threads;
            int start = thread_id * points_per_thread;
            int end = thread_id * points_per_thread + points_per_thread - 1;

            vector<Point> slice_x = slice_vec(local_x, start, end);
            vector<Point> slice_y = slice_x;
            sort(slice_y.begin(), slice_y.end(), [](Point const & a, Point const & b) -> bool
            { return a.y < b.y; } 
            );

            delta[thread_id] = min_dist(slice_x, slice_y); // each thread calculates their local delta
            local_delta = delta[thread_id];
            #pragma omp barrier

            if (thread_id % 2 == 0) {
                double delta_thread = min(delta[thread_id], delta[thread_id + 1]);
                double middle = (slice_x[slice_x.size() - 1].x + local_x[end + 1].x) / 2.0; 
                vector<Point> strip = compute_strip(slice_y, delta_thread, middle);
                strips[thread_id] = strip;

            } else {
                double delta_thread = min(delta[thread_id], delta[thread_id - 1]);
                double middle = (slice_x[0].x + local_x[start - 1].x) / 2.0; 
                vector<Point> strip = compute_strip(slice_y, delta_thread, middle);
                strips[thread_id] = strip;
            }
            #pragma omp barrier

            if (thread_id % 2 == 0) {
                vector<Point> strip = merge_two_vectors(strips[thread_id], strips[thread_id + 1], false);
                double tmp = DBL_MAX;

                // for all points in strip check the next 6 points down the strip
                for (int i = 0; i < strip.size(); i++) 
                    for (int j = i + 1; j < i + 7 && j < strip.size(); j++) 
                        tmp = min(tmp, sqrt(dist(strip[i], strip[j])));

                local_delta = min(tmp, local_delta);
            }
        }
    }

    if (rank % 2) {   // right side
        MPI_Send(&local_delta, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);   // send local delta
        MPI_Recv(&local_delta, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // receive min delta
        
        double last_left; 
        MPI_Recv(&last_left, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // receive last point on the left hand side
        double middle_point = (local_x[0].x - last_left) / 2.0;
        MPI_Send(&middle_point, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);   // send middle point x-coordinate

        vector<Point> strip = compute_strip(local_y, local_delta, middle_point);
        int strip_size = strip.size();
        MPI_Send(&strip_size, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD);   // send strip size
        MPI_Send(strip.data(), strip.size(), dt_point, rank - 1, 0, MPI_COMM_WORLD);   // send strip points 

        double infinite = DBL_MAX;
        MPI_Reduce(&infinite, &infinite, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    } else {    // left side
        double ext_delta;
        MPI_Recv(&ext_delta, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // receive external delta
        local_delta = min(local_delta, ext_delta);
        MPI_Send(&local_delta, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);   // send local delta

        MPI_Send(&local_x[local_x.size() - 1].x, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);   // send last point on this side
        double middle_point;
        MPI_Recv(&middle_point, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // receive middle point x-coordinate      

        vector<Point> strip;
        {
            vector<Point> strip_one = compute_strip(local_y, local_delta, middle_point);
            int msg_size;
            MPI_Recv(&msg_size, 1 , MPI_INT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // receive the second strips size    
            vector<Point> strip_two(msg_size);
            MPI_Recv(strip_two.data(), msg_size, dt_point, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // receive the second strip      
            strip = merge_two_vectors(strip_one, strip_two, false);
        }

        // for all points in strip check the next 6 points down the strip
        for (int i = 0; i < strip.size(); i++) 
            for (int j = i + 1; j < i + 7 && j < strip.size(); j++) 
                local_delta = min(local_delta, sqrt(dist(strip[i], strip[j])));

        double result;
        MPI_Reduce(&local_delta, &result, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        if (rank == 0)
            cout << "RESULT = " << result << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

/**
 * @brief O(n log n) divide-and-conquer algorithm to compute the distance between any two points.
 * 
 * @param x_sorted The set of all points sorted by the x-coordinate 
 * @param y_sorted The set of all points sorted by the y-coordinate 
 * 
 * @return The minimal distance between any two points
*/
double min_dist(vector<Point>& x_sorted, vector<Point>& y_sorted) {
    int n = x_sorted.size();
    if (n < 10) 
        return brute_force(x_sorted);
    
    Point median = x_sorted[n / 2];
    vector<Point> x_sorted_left = vector<Point>(x_sorted.begin(), x_sorted.begin() + n / 2);
    vector<Point> x_sorted_right = vector<Point>(x_sorted.begin() + n / 2 + 1, x_sorted.end());
    vector<Point> y_sorted_left;
    vector<Point> y_sorted_right;

    // assign points left of the line to ysorted_left and the others to y_sorted_right
    for (auto p : y_sorted)
        if (p.x <= median.x)
            y_sorted_left.push_back(p);
        else   
            y_sorted_right.push_back(p);

    // recursively solve left and right     
    double delta_1 = min_dist(x_sorted_left, y_sorted_left); 
    double delta_2 = min_dist(x_sorted_right, y_sorted_right); 
    double delta = min(delta_1, delta_2);

    // put points whose distance is at most delta from the vertical strip inside strip
    vector<Point> strip;
    for (auto p : y_sorted) 
        if (p.x < median.x + delta && p.x > median.x - delta)
            strip.push_back(p);

    // for all points in strip check the next 6 points down the strip
    for (int i = 0; i < strip.size(); i++) 
        for (int j = i + 1; j < i + 7 && j < strip.size(); j++) 
            delta = min(delta, sqrt(dist(strip[i], strip[j])));
        
    return delta;
}

/**
 * @brief Compute points within delta of middle from a set of points
 * 
 * @param y_sorted The set of points
 * @param delta the size of the strip left and right
 * @param middle the x-coordinate of the center of the strip
 * 
 * @return Vector of points in the strip
*/
vector<Point> compute_strip(vector<Point>& y_sorted, double delta, double middle) {
    vector<Point> strip;
        for (auto p : y_sorted) 
            if (fabs(middle - p.x) < delta)
                strip.push_back(p);

    return strip;
}

/**
 * @brief sort a vector using 1 process and many threads
 * 
 * @param A The vector to sort
 * @param threads The number of threads
*/
void sort_using_threads(vector<Point> &A, int threads) {
    #pragma omp parallel num_threads(threads)
    {
        int thread_id = omp_get_thread_num();
        int n_points_to_sort = A.size() / threads;

        if (thread_id != threads - 1)
            sort(A.begin() + thread_id * n_points_to_sort, A.begin() + thread_id * n_points_to_sort + n_points_to_sort, 
            [](Point const & a, Point const & b) -> bool
            { return a.y < b.y; } 
            );
        else 
            sort(A.begin() + thread_id * n_points_to_sort, A.end(), 
            [](Point const & a, Point const & b) -> bool
            { return a.y < b.y; } 
            );
    }

    vector<Point> A_copy(A.size());

    #pragma omp parallel num_threads(threads)
    {
        int thread_id = omp_get_thread_num();
        bool turn = true; // used to decide whether to write in A or A_copy
        int msg_size = A.size() / threads, t_size = threads;

        while (t_size != 1) {
            if (thread_id < t_size && thread_id < t_size / 2) {
                if (!turn) {
                    merge_two_vectors_in_place(A_copy, A, thread_id * msg_size, (thread_id + t_size / 2) * msg_size, thread_id * msg_size * 2, msg_size, false);
                } else { 
                    merge_two_vectors_in_place(A, A_copy, thread_id * msg_size, (thread_id + t_size / 2) * msg_size, thread_id * msg_size * 2, msg_size, false);
                }
            }
            t_size /= 2;
            msg_size *= 2;
            turn = !turn;
            #pragma omp barrier 
        }
        if (!turn && thread_id == 0)
            A = A_copy;
    }
}

/**
 * @brief Sort a vector of points in parallel using both MPI and OPENMP by the x-coordinate.
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

    vector<Point> local_points_copy(local_points.size());
    #pragma omp parallel num_threads(threads)
    {
        int thread_id = omp_get_thread_num();
        bool turn = true; // used to decide whether to write in local_points or local_points_copy
        int msg_size = local_points.size() / threads, t_size = threads;

        while (t_size != 1) {
            if (thread_id < t_size && thread_id < t_size / 2) {
                if (!turn) {
                    merge_two_vectors_in_place(local_points_copy, local_points, thread_id * msg_size, (thread_id + t_size / 2) * msg_size, thread_id * msg_size * 2, msg_size, true);
                } else { 
                    merge_two_vectors_in_place(local_points, local_points_copy, thread_id * msg_size, (thread_id + t_size / 2) * msg_size, thread_id * msg_size * 2, msg_size, true);
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
void merge_two_vectors_in_place(vector<Point> &src, vector<Point> &dst, int i, int j, int k, int n, bool sort_by_x) {
    int i_max = i + n, j_max = j + n;

    while (i < i_max && j < j_max) {
        if (sort_by_x) {
            if (src[i].x < src[j].x)
                dst[k++] = src[i++];
            else
                dst[k++] = src[j++];
        } else {
            if (src[i].y < src[j].y)
                dst[k++] = src[i++];
            else
                dst[k++] = src[j++];
        }
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
    int m = B.size();
    vector<Point> C(n + m);
    int i = 0, j = 0, k = 0;

    while (i < n && j < m) {
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
        while (j < m)
            C[k++] = B[j++];
    else 
        while (i < n)
            C[k++] = A[i++];

    return C;
}

/**
 * @brief Take a slice of A from start to end
 * 
 * @param A The vector to take the slice from
 * @param start starting index
 * @param end ending index
*/
vector<Point> slice_vec(vector<Point> &A, int start, int end) {
    vector<Point> B(end - start + 1);
    int k = 0;
    for (int i = start; i <= end; i++) 
        B[k++] = A[i];
    return B;
}

/**
 * @brief Distance between two points in 2 dimensions.
 * 
 * @param p1 The first point
 * @param p2 The second point
 * 
 * @return The distance between p1 and p2
 * @note In order to make the computation lighter, the square root is not executed. Thus, if this value needs to be actually used, it needs to be 'squared rooted'.
*/
inline double dist(Point& p1, Point& p2) {
    double delta_x = fabs(p1.x - p2.x);
    double delta_y = fabs(p1.y - p2.y);
    return delta_x * delta_x + delta_y * delta_y;
}

/**
 * @brief computes the minimal distance between any two points
 * 
 * @param points The set of points to check
 * 
 * @return The minimal distance between any two points (without square root)
*/
double brute_force(vector<Point>& points) {
    double min_dist = dist(points[0], points[1]);
    int n = points.size();

    for (int i = 0; i < n; i++) 
        for (int j = i + 1; j < n; j++) 
            min_dist = min(min_dist, dist(points[i], points[j]));
        
    return sqrt(min_dist); 
}