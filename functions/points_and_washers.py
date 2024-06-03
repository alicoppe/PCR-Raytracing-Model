import numpy as np
from functions.helper_functions import pcr_cone_height, x_value_at

def pcr_tube_to_washers(n, v):
    """
    Converts PCR tube volume and number of washers into radius values from top to bottom.

    Parameters:
    n: Number of washers
    v: Volume of PCR tube

    Returns:
    r: Array of radius values
    dz: Incremental heights of each washer
    height: Height of the volume with respect to the z-axis
    """
    height = pcr_cone_height(v)
    m, b = 6.8249, -7.0860
    a, c = 1.3080, 1.120
    elevation_of_water = 0.721
    
    diff = height - elevation_of_water
    dz = diff / n
    r = np.zeros(n)
    
    for i in range(n):
        spacing = dz * (i + 0.5)
        z = height - spacing
        if z >= (c + elevation_of_water):
            r[i] = x_value_at(z, m, b)
        else:
            z2 = z - (c + elevation_of_water)
            r[i] = np.sqrt((1 - (z2 ** 2 / c ** 2)) * a ** 2)
    
    return r, dz, height

def excitation_power_array2(r, dz, extinction_coeff, c, p_i):
    """
    Calculates the excitation power at each washer level.

    Parameters:
    r: Radius array
    dz: Incremental heights
    extinction_coeff: Extinction coefficient
    c: Constant
    p_i: Initial power

    Returns:
    excitation_powers: Array of excitation powers
    """
    excitation_powers = np.zeros(len(r))
    for i in range(len(r)):
        P_d = p_i * np.exp(-extinction_coeff * (c * dz))
        P_abs = p_i - P_d
        excitation_powers[i] = P_abs
        p_i = P_d
    
    return excitation_powers

def volume_of_washers(r, dz):
    """
    Calculates the total volume of washers.

    Parameters:
    r: Radius array
    dz: Incremental heights

    Returns:
    volume: Total volume
    """
    volume = np.sum(np.pi * r ** 2 * dz)
    return volume

def point_distribution(r, dz, total_points):
    """
    Distributes points based on the volume of washers.

    Parameters:
    r: Radius array
    dz: Incremental heights
    total_points: Total number of points to distribute

    Returns:
    num_points: Array of point distribution per washer
    """
    volume = volume_of_washers(r, dz)
    remaining_points = total_points
    remaining_volume = volume
    num_points = np.zeros(len(r), dtype=int)
    
    for i in range(len(r)):
        this_vol = np.pi * r[i] ** 2 * dz
        percent_of_remaining_volume = this_vol / remaining_volume
        num_points[i] = round(percent_of_remaining_volume * remaining_points)
        remaining_points -= num_points[i]
        remaining_volume -= this_vol
    
    return num_points

def random_points_circle(n, r):
    """
    Generates random points within a circle.

    Parameters:
    n: Number of points
    r: Radius of the circle

    Returns:
    x, y: Arrays of x and y coordinates of points
    """
    radius = np.sqrt(np.random.rand(n)) * r
    thetas = np.random.rand(n) * 2 * np.pi
    x = radius * np.cos(thetas)
    y = radius * np.sin(thetas)
    return x, y

def generate_washer_points(r, dz, height, total_points):
    """
    Generates points within washers.

    Parameters:
    r: Radius array
    dz: Incremental heights
    height: Height of the volume
    total_points: Total number of points to generate

    Returns:
    points: Array of generated points with (x, y, z) coordinates
    """
    points = np.zeros((total_points, 3))
    num_points = point_distribution(r, dz, total_points)
    counter = 0
    
    for i in range(len(r)):
        z = height - dz * (i + 0.5)
        for j in range(num_points[i]):
            x, y = random_points_circle(1, r[i])
            points[counter,0] = x
            points[counter,1] = y
            points[counter,2] = z
            counter += 1
    
    return points
