#include <bits/stdc++.h>
#include <sys/time.h>

using namespace std;

class Point {
public:
    double x; 
    double y;
};

vector<Point> generate_points(int);
double closest_pair(vector<Point>&);
double dist(Point &, Point&);

int main() {
    cout << "INPUT SIZE \t TIME\n";

    for (int i = 10; i < 19; i++) {
        vector<Point> points = generate_points(1 << i);
        struct timeval start, end;

        gettimeofday(&start, NULL);
        double min_dist = closest_pair(points);
        gettimeofday(&end, NULL);

        long long int time_usec = ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));
        double time_sec = (double)time_usec / 1000000.0;
        cout << (1 << i) << "\t\t" << time_sec << "\n";
    }
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

double closest_pair(vector<Point>& points) {
    double min_dist = dist(points[0], points[1]);
    int n = points.size();

    for (int i = 0; i < n; i++) 
        for (int j = i + 1; j < n; j++) 
            min_dist = min(min_dist, dist(points[i], points[j]));
        
    return sqrt(min_dist); 
}

inline double dist(Point& p1, Point& p2) {
    double delta_x = fabs(p1.x - p2.x);
    double delta_y = fabs(p1.y - p2.y);
    return delta_x * delta_x + delta_y * delta_y;
}

