#include <bits/stdc++.h>
using namespace std;

class Point {
public:
    double x; 
    double y;
};

int N;

vector<Point> generate_points(int); 
double brute_force(vector<Point>& points);
inline double dist(Point& p1, Point& p2);
double min_dist(vector<Point>& x_sorted, vector<Point>& y_sorted);

int main(int argc, char** argv) {
    if (argc == 1)
        return 0;
    N = atoi(argv[1]);

    vector<Point> x_sorted = generate_points(N);
    vector<Point> y_sorted = x_sorted;
    
    sort(x_sorted.begin(), x_sorted.end(), [](Point const & a, Point const & b) -> bool
        { return a.x < b.x; } 
    );
    sort(y_sorted.begin(), y_sorted.end(), [](Point const & a, Point const & b) -> bool
        { return a.y < b.y; } 
    );

    cout << "MIN DIST = " << min_dist(x_sorted, y_sorted) << endl;
}

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
    for (auto &p : y_sorted)
        if (p.x < median.x)
            y_sorted_left.push_back(p);
        else
            y_sorted_right.push_back(p);

    // recursively solve left and right     
    double delta_1 = min_dist(x_sorted_left, y_sorted_left); 
    double delta_2 = min_dist(x_sorted_right, y_sorted_right); 
    double delta = min(delta_1, delta_2);
    
    // put points whose distance is at most delta from the vertical strip inside strip
    vector<Point> strip;
    for (auto &p : y_sorted)   
        if (p.x < median.x + delta && p.x > median.x - delta)
            strip.push_back(p);
    
    // for all points in strip check the next 6 points down the strip
    for (int i = 0; i < strip.size(); i++) 
        for (int j = i + 1; j < i + 7 && j < strip.size(); j++) 
            delta = min(delta, dist(strip[i], strip[j]));
        
    return delta;
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