#include <bits/stdc++.h>
using namespace std;

class Point {
public:
    double x; 
    double y;
};

vector<Point> generate_points(int n);
inline double dist(Point& p1, Point& p2);
double brute_force(vector<Point>& points);
ifstream in("in.txt");

int main() {
    vector<Point> points(1024);

    for (int i = 0; i < 1024; i++) {
        Point p;
        in >> p.x >> p.y;
        points[i] = p;
    }

    cout << "RESULT = " << brute_force(points) << endl;
}

double brute_force(vector<Point>& points) {
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