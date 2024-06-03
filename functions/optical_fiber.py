import numpy as np

def cone_acceptance_angle(acceptance_angle, mid, radius):
    """
    Generates points for a cone based on the acceptance angle.

    Parameters:
    acceptance_angle: Angle of acceptance in radians
    mid: Midpoint coordinates as a list or array [x, y, z]
    radius: Radius of the base of the cone

    Returns:
    x, y, z: Arrays of x, y, z coordinates of the points on the cone
    """
    rh_4 = acceptance_angle
    th_4 = np.linspace(0, 2 * np.pi, int(2 * np.pi / 0.05) + 1)
    rad_4 = np.linspace(0, 10, int(10 / 0.3) + 1)
    x_start = (radius / np.tan(acceptance_angle)) + mid[0]
    THETA_4, RAD_4 = np.meshgrid(th_4, rad_4)
    x = -RAD_4 * np.cos(rh_4) + x_start
    y = RAD_4 * np.sin(rh_4) * np.cos(THETA_4) + mid[1]
    z = RAD_4 * np.sin(THETA_4) * np.sin(rh_4) + mid[2]
    return x, y, z

def optical_fiber_head(mid, r):
    """
    Generates points for the head of an optical fiber.

    Parameters:
    mid: Midpoint coordinates as a list or array [x, y, z]
    r: Radius of the optical fiber head

    Returns:
    x, y, z: Arrays of x, y, z coordinates of the points on the optical fiber head
    """
    theta = np.linspace(0, 2 * np.pi, int(2 * np.pi / 0.05) + 1)
    x1 = np.linspace(0, 0.5, int(0.5 / 0.1) + 1)
    THETA, X = np.meshgrid(theta, x1)
    x = X + mid[0]
    y = r * np.cos(THETA) + mid[1]
    z = r * np.sin(THETA) + mid[2]
    return x, y, z
