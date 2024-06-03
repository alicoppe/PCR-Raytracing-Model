import numpy as np
from functions.helper_functions import lin_eq, distance_2d

def ellipsoid1(a, b, c, theta, rho, x_i, y_i, z_i):
    """
    Generates coordinates for an ellipsoid.

    Parameters:
    a, b, c: Ellipsoid parameters
    theta: Spherical coordinate bounds (array of two elements)
    rho: Spherical coordinate bounds (array of two elements)
    x_i, y_i, z_i: Initial coordinates

    Returns:
    x, y, z: Coordinates of the ellipsoid
    """
    th = np.linspace(theta[0], theta[1], 11)
    rh = np.linspace(rho[0], rho[1], 11)
    THETA, RHO = np.meshgrid(th, rh)
    x = a * np.sin(RHO) * np.cos(THETA) + x_i
    y = b * np.sin(THETA) * np.sin(RHO) + y_i
    z = c * np.cos(RHO) + z_i
    return x, y, z

def vertical_cone(rad_1, rad_2, height, x_i, y_i, z_i):
    """
    Generates coordinates for a vertical cone.

    Parameters:
    rad_1: Bottom radius
    rad_2: Top radius
    height: Distance from rad_1 to rad_2
    x_i, y_i, z_i: Initial coordinates

    Returns:
    x, y, z: Coordinates of the cone
    """
    m, b = lin_eq([rad_1, z_i], [rad_2, z_i + height])
    rho = np.arctan(1 / m)
    start = distance_2d(0, rad_1, b, z_i)
    finish = distance_2d(0, rad_2, b, z_i + height)
    theta = np.linspace(0, 2 * np.pi, 41)
    radius = np.linspace(start, finish, 31)
    THETA, RADIUS = np.meshgrid(theta, radius)
    x = RADIUS * np.sin(rho) * np.cos(THETA) + x_i
    y = RADIUS * np.sin(THETA) * np.sin(rho) + y_i
    z = RADIUS * np.cos(rho) + b
    return x, y, z

def vertical_cylinder(radius, height, x_i, y_i, z_i):
    """
    Generates coordinates for a vertical cylinder.

    Parameters:
    radius: Radius of the cylinder
    height: Height of the cylinder
    x_i, y_i, z_i: Initial coordinates

    Returns:
    x, y, z: Coordinates of the cylinder
    """
    theta = np.linspace(0, 2 * np.pi, 41)
    heights = np.linspace(z_i, height + z_i, 31)
    THETA, HEIGHTS = np.meshgrid(theta, heights)
    x = radius * np.cos(THETA) + x_i
    y = radius * np.sin(THETA) + y_i
    z = HEIGHTS
    return x, y, z
