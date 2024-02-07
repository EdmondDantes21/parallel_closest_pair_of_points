#include <bits/stdc++.h>
#include <mpi.h>
#include <sys/time.h>

using namespace std;

class Point {
public:
    double x; 
    double y;
};

vector<Point> generate_points(int); 
double min_dist(vector<Point>& x_sorted, vector<Point>& y_sorted);
double brute_force(vector<Point>& points);
inline double dist(Point& p1, Point& p2);
vector<Point> compute_strip(vector<Point>& y_sorted, double delta, double middle);
vector<Point> merge_two_vectors(vector<Point> &A, vector<Point> &B, bool sort_by_x);
void parallel_merge_sort(vector<Point> &points, int world_rank, int world_size, MPI_Datatype dt_point);

int N;

int main(int argc, char **argv) {
    if (argc == 1)
        return 0;
    N = atoi(argv[1]);
    struct timeval start,end;

    gettimeofday(&start, NULL);

    MPI_Init(NULL, NULL);
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // custom data type
    MPI_Datatype dt_point;
    MPI_Type_contiguous(2, MPI_DOUBLE, &dt_point);
    MPI_Type_commit(&dt_point);

    // root process holds the data
    vector<Point> points;
    if (world_rank == 0)
        points = generate_points(N);
    
    parallel_merge_sort(points, world_rank, world_size, dt_point);

    vector<Point> local_x(N / world_size);
    MPI_Scatter(points.data(), N / world_size, dt_point, local_x.data(), N / world_size, dt_point, 0, MPI_COMM_WORLD);
    vector<Point> local_y = local_x;
    sort(local_y.begin(), local_y.end(), [](Point const & a, Point const & b) -> bool
        { return a.x < b.x; } 
    );

    if (world_rank == 0)
        points.clear();

    double local_delta = min_dist(local_x, local_y); // each process calculates their local delta

    if (world_rank % 2) {   // right side
        MPI_Send(&local_delta, 1, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD);   // send local delta
        MPI_Recv(&local_delta, 1, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // receive min delta
        
        double last_left; 
        MPI_Recv(&last_left, 1, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // receive last point on the left hand side
        double middle_point = (local_x[0].x - last_left) / 2.0;
        MPI_Send(&middle_point, 1, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD);   // send middle point x-coordinate

        vector<Point> strip = compute_strip(local_y, local_delta, middle_point);
        int strip_size = strip.size();
        MPI_Send(&strip_size, 1, MPI_INT, world_rank - 1, 0, MPI_COMM_WORLD);   // send strip size
        MPI_Send(strip.data(), strip.size(), dt_point, world_rank - 1, 0, MPI_COMM_WORLD);   // send strip points 

        double infinite = DBL_MAX;
        MPI_Reduce(&infinite, &infinite, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    } else {    // left side
        double ext_delta;
        MPI_Recv(&ext_delta, 1, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // receive external delta
        local_delta = min(local_delta, ext_delta);
        MPI_Send(&local_delta, 1, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD);   // send local delta

        MPI_Send(&local_x[local_x.size() - 1].x, 1, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD);   // send last point on this side
        double middle_point;
        MPI_Recv(&middle_point, 1, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // receive middle point x-coordinate      

        vector<Point> strip;
        {
            vector<Point> strip_one = compute_strip(local_y, local_delta, middle_point);
            vector<Point> strip_two;
            int msg_size;
            MPI_Recv(&msg_size, 1 , MPI_INT, world_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // receive the second strips size    
            MPI_Recv(strip_two.data(), msg_size, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // receive the second strip      
            strip = merge_two_vectors(strip_one, strip_two, false);
        }

        // for all points in strip check the next 6 points down the strip
        for (int i = 0; i < strip.size(); i++) 
            for (int j = i + 1; j < i + 7 && j < strip.size(); j++) 
                local_delta = min(local_delta, sqrt(dist(strip[i], strip[j])));

        double result;
        MPI_Reduce(&local_delta, &result, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        gettimeofday(&end, NULL);

        if (world_rank == 0) {
            gettimeofday(&end, NULL);
            long long int time_usec = ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));
            double time_sec = (double)time_usec / 1000000.0;
            cout << "N = " << N << endl;
            cout << "T = " << time_sec << endl;
        }  
    }
    
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
void parallel_merge_sort(vector<Point> &points, int rank, int size, MPI_Datatype dt_point) {
    vector<Point> local_points(N / size);
    MPI_Scatter(points.data(), N / size, dt_point, local_points.data(), N / size, dt_point, 0, MPI_COMM_WORLD);
    sort(local_points.begin(), local_points.end(), [](Point const & a, Point const & b) -> bool
        { return a.x < b.x; } 
    );

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
        for (auto &p : y_sorted) 
            if (p.x - delta > middle && p.x + delta < middle)
                strip.push_back(p);
    return strip;
}
/**
 * @brief Merge vectors A and B while keeping the y-coordinates order
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