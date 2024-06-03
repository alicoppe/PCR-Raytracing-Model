import numpy as np

def find_angle(mid, point, thickness, d):
    dist = distance_1d(point, mid)
    half_h = thickness / 2
    if dist < half_h:
        h_1 = half_h + dist
        h_2 = half_h - dist
        angle = np.arctan(h_1 / d) + np.arctan(h_2 / d)
    elif dist > half_h:
        angle = np.arctan((half_h + dist) / d) - np.arctan((dist - half_h) / d)
    else:
        angle = np.arctan(thickness / d)
    return angle

def distance_1d(start, finish):
    return abs(finish - start)

def distance_2d(xi, xf, yi, yf):
    return np.sqrt((xf - xi) ** 2 + (yf - yi) ** 2)

def distance_3d(p1, p2):
    return np.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2 + (p2[2] - p1[2]) ** 2)

def lin_eq(p1, p2):
    m = (p2[1] - p1[1]) / (p2[0] - p1[0])
    b = p1[1] - m * p1[0]
    return m, b

def integral_of_circle(xi, xf, r):
    return (r ** 2) * 0.5 * (np.arcsin(xf / r) - np.arcsin(xi / r) + 0.5 * (np.sin(2 * np.arcsin(xf / r)) - np.sin(2 * np.arcsin(xi / r))))

def integral_of_line(xi, xf, m, b):
    return (m / 2) * (xf ** 2 - xi ** 2) + b * (xf - xi)

def y_value_at(x, m, b):
    return m * x + b

def x_value_at(y, m, b):
    return (y - b) / m

def line_intersection(m1, b1, m2, b2):
    return (b2 - b1) / (m1 - m2)

def pcr_cone_height(v):
    a1 = 1.3080
    b1 = 1.3080
    c1 = 1.120
    m = 6.8249
    b = -7.0860
    h1 = 0.721 + 1.120

    k1 = np.pi / m
    k2 = (2 * np.pi * b) / m

    v = v - (4 / 3) * np.pi * a1 * b1 * c1
    c = k2 * h1 - k1 * (h1 ** 2) - v
    p = [k1, -k2, c]
    q = np.roots(p)
    height = q[q > 0]
    return height[0] if height.size > 0 else None

def intersection_line_with_circle(m, b, r):
    a = m ** 2 + 1
    b1 = 2 * b * m
    c = b ** 2 - r ** 2

    x = (-b1 + np.sqrt(b1 ** 2 - 4 * a * c)) / (2 * a)
    return x


def intersection_line_with_circle2(m, b, r):
    a = m ** 2 + 1
    b1 = 2 * b * m
    c = b ** 2 - r ** 2

    roots = np.roots([a, b1, c])
    
    
    if np.isnan(roots).any() or np.isinf(roots).any():
        # Debug statements
        print(f"intersection_line_with_circle2: m: {m}, b: {b}, r: {r}, roots: {roots}")
        raise ValueError("Roots contain NaN or infinite values")
    
    return np.max(roots)

def find_solid_angle(array_point, array_midpoint, h, w):
    d = abs(array_midpoint[0] - array_point[0])
    theta_max = find_angle(array_point[2], array_midpoint[2], h, d)
    phi_max = find_angle(array_point[1], array_midpoint[1], w, d)
    omega = phi_max * (-np.cos(theta_max) + 1)
    return omega
