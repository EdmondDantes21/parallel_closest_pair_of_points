#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <vector>
#include <math.h>
#include <fstream>
#include <sys/time.h>
#include <cstdlib>

using namespace std;

class Point {
public:
    double x; 
    double y;
};

vector<Point> generate_points(int n);
double closest_pair(vector<Point>&);
double dist(Point &, Point&);

int main(int argc, char** argv)
{
    if (argc == 0)
        return 0;

    int N = atoi(argv[1]);
    struct timeval start,end;

    vector<Point> points = generate_points(N);

    MPI_Init(NULL, NULL);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0)
        gettimeofday(&start, NULL);

    double min_dist = dist(points[0], points[1]);
    for (int i = world_rank; i < N; i += world_size) 
        for (int j = i + 1; j < N; j++) 
            min_dist = min(min_dist, dist(points[i], points[j]));
        
    min_dist = sqrt(min_dist);

    double result;
    MPI_Reduce(&min_dist, &result, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

    if (world_rank == 0) {
        gettimeofday(&end, NULL);
        long long int time_usec = ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));
        double time_sec = (double)time_usec / 1000000.0;
        cout << "T = " << time_sec << endl;
    }

    MPI_Finalize();
    return 0;
}

vector<Point> generate_points(int n) {
    vector<Point> points(n);

    srandom(time(NULL));
    const long max_rand = 1000000L;
    double lower_bound = -10000.0;
    double upper_bound = 10000.0;

    for (int i = 0; i < n; i++) {
        Point p;
        p.x = lower_bound + (upper_bound - lower_bound) * (random() % max_rand) / max_rand;
        p.y = lower_bound + (upper_bound - lower_bound) * (random() % max_rand) / max_rand;
        points[i] = p;
    }
    return points;
}

inline double dist(Point& p1, Point& p2) {
    double delta_x = fabs(p1.x - p2.x);
    double delta_y = fabs(p1.y - p2.y);
    return delta_x * delta_x + delta_y * delta_y;
}