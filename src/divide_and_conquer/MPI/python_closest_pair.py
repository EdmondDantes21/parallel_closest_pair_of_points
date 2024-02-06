import math

def dist(p1, p2):
    return math.sqrt(((p2[1]-p1[1])**2)+((p2[0]-p1[0])**2))

def closest_brute_force(points):
    min_dist = float("inf")
    p1 = None
    p2 = None
    for i in range(len(points)):
        for j in range(i+1, len(points)):
            d = dist(points[i], points[j])
            if d < min_dist:
                min_dist = d
                p1 = points[i]
                p2 = points[j]
    return p1, p2, min_dist


def rec(xsorted, ysorted):
    n = len(xsorted)
    if n <= 3:
        return closest_brute_force(xsorted)
    else:
        midpoint = xsorted[n//2]    # median point
        xsorted_left = xsorted[:n//2]   # points on the left sorted by x
        xsorted_right = xsorted[n//2:]  # points on the right sorted by x
        ysorted_left = []   # points on the left sorted by y
        ysorted_right = []  # points on the right sorted by y

        # assign points left of the line to ysorted_left and the others to y_sorted_right
        for point in ysorted:
            ysorted_left.append(point) if (point[0] <= midpoint[0]) else ysorted_right.append(point)
        
        # solve the problem recursively for the left and right side
        (p1_left, p2_left, delta_left) = rec(xsorted_left, ysorted_left)
        (p1_right, p2_right, delta_right) = rec(xsorted_right, ysorted_right)

        # find the minimum delta between the two sides
        (p1, p2, delta) = (p1_left, p2_left, delta_left) if (delta_left < delta_right) else (p1_right, p2_right, delta_right)

        # in_band contains points in the famous band, sorted by y
        in_band = [point for point in ysorted if midpoint[0]-delta < point[0] < midpoint[0] + delta]

        for i in range(len(in_band)):
            for j in range(i+1, min(i+7, len(in_band))): # at most 6 points down
                d = dist(in_band[i], in_band[j])
                if d < delta:
                    (p1, p2, delta) = (in_band[i], in_band[j], d)
                    
        return p1, p2, delta


def closest(points):
    xsorted = sorted(points, key=lambda point: point[0]) # all points sorted by x
    ysorted = sorted(points, key=lambda point: point[1]) # all points sorted by y
    return rec(xsorted, ysorted)