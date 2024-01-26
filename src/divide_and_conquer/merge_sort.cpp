#include <bits/stdc++.h>
#include <iostream>
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

int main() {
    vector<Point> points = generate_points(10);
    merge_sort(points, 0, 9);
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