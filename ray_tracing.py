import numpy as np
from functions.points_and_washers import pcr_tube_to_washers, generate_washer_points, excitation_power_array2, point_distribution
from functions.helper_functions import distance_2d, find_solid_angle, distance_3d, lin_eq, y_value_at, intersection_line_with_circle2, line_intersection


def find_total_power(total_points, num_vectors, volume, mid, radius, acceptance_angle, excitation_power, quantum_yield, c, extinction_coeff):
    # Parameters for parametric modelling of the volumes
    rho1 = 0.1455
    rho2 = 0.1471
    p1 = [1.3080, 1.8410]  # Parameters for the first ellipsoid
    p2 = [1.4695, 1.270]   # Parameters for the second ellipsoid

    h = (radius * 2) * 1.5  # Height of the region to be modeled
    w = h  # Width is equal to height
    num_washers = round(total_points / 60)  # Number of washers to generate
    r, dz, height = pcr_tube_to_washers(num_washers, volume)  # Generate washer dimensions
    points = generate_washer_points(r, dz, height, total_points)  # Generate points within washers
    power_factors = np.zeros(num_washers)  # Initialize power factors array
    excitation_powers = excitation_power_array2(r, dz / 10, extinction_coeff, c, excitation_power)  # Calculate excitation powers
    num_points = point_distribution(r, dz, total_points)  # Distribute points across washers

    rho = acceptance_angle * 1.3  # Adjusted acceptance angle
    x_start = (radius / np.tan(rho)) + mid[0]  # Starting x-coordinate for point calculation

    counter = 0
    for i in range(num_washers):
        power_factors_washer = np.zeros(num_points[i])  # Initialize power factors for each washer
        for j in range(num_points[i]):
            x = -(points[counter, 0] - x_start)  # Calculate x distance from start
            dist = x * np.tan(rho)  # Calculate radial distance based on acceptance angle
            if distance_2d(mid[1], points[counter, 1], mid[2], points[counter, 2]) <= dist:  # Check if point is within acceptance cone
                initial_vectors = generate_vectors(h, w, mid, points[counter, :], num_vectors)  # Generate initial vectors
                initial_startpoints = np.ones((num_vectors, 3)) * points[counter, :]  # Set initial start points

                # Calculate intersections and refractions through ellipsoids
                initial_endpoints, normal_vectors = vect_cone_intersection(initial_startpoints, initial_vectors, rho1, p1)
                refracted_1 = refracted_vectors(initial_vectors, normal_vectors, 1.333, 1.49)
                initial_endpoints_2, normal_vectors_2 = vect_cone_intersection(initial_endpoints, refracted_1, rho2, p2)
                refracted_2 = refracted_vectors(refracted_1, normal_vectors_2, 1.49, 1)
                endpoints = vector_endpoints(refracted_2, initial_endpoints_2, mid[0])

                accepted_indices = optical_fiber_acceptance_quick(radius, acceptance_angle, mid, refracted_2, endpoints)  # Find accepted vectors
                a_startpoints = initial_endpoints[accepted_indices > 0, :]  # Filter accepted start points

                if a_startpoints.size > 0:
                    solid_angle = find_solid_angle(points[counter, :], mid, h, w)  # Calculate solid angle
                    if extinction_coeff > 0:
                        # Calculate power factors for the point with extinction coefficient
                        power_factors_washer[j] = point_average_power2(excitation_powers[i], points[counter, :], num_vectors, solid_angle, a_startpoints, quantum_yield, c, extinction_coeff)
                    else:
                        # Calculate power factors for the point without extinction coefficient
                        vector_ratio = len(a_startpoints) / num_vectors
                        power_factors_washer[j] = point_average_power3(excitation_powers[i], vector_ratio, solid_angle, quantum_yield)
            
            counter += 1  # Increment point counter
        power_factors[i] = np.mean(power_factors_washer)  # Average power factors for the washer

    ratio_total_power = np.sum(power_factors)  # Sum of power factors to get total power ratio
    return ratio_total_power

def vector_plot(origin, num_vectors, mid, radius, acceptance_angle):
    # Parameters for parametric modelling of the volumes
    rho1 = 0.1455
    rho2 = 0.1471
    p1 = [1.3080, 1.8410]  # Parameters for the first ellipsoid
    p2 = [1.4695, 1.270]   # Parameters for the second ellipsoid

    height = (radius * 2) * 1.5  # Height of the region to be modeled
    width = height  # Width is equal to height

    initial_vectors = generate_vectors(height, width, mid, origin, num_vectors)  # Generate initial vectors
    initial_startpoints = np.ones((num_vectors, 3)) * origin  # Set initial start points
    initial_endpoints, normal_vectors = vect_cone_intersection(initial_startpoints, initial_vectors, rho1, p1)  # Calculate intersections with first ellipsoid

    initial = np.hstack((initial_startpoints, initial_endpoints - initial_startpoints))  # Combine start and end points for initial vectors

    refracted_1 = refracted_vectors(initial_vectors, normal_vectors, 1.33, 1.49)  # Refract vectors at first interface
    refracted_endpoints_1, normal_vectors_2 = vect_cone_intersection(initial_endpoints, refracted_1, rho2, p2)  # Calculate intersections with second ellipsoid
    refracted_2 = refracted_vectors(refracted_1, normal_vectors_2, 1.49, 1)  # Refract vectors at second interface
    refracted_endpoints_2 = vector_endpoints(refracted_2, refracted_endpoints_1, mid[0])  # Calculate final endpoints
    initial_2 = np.hstack((initial_endpoints, refracted_endpoints_1 - initial_endpoints))  # Combine start and end points for refracted vectors

    # Determine accepted and not accepted vectors
    a_end, na_end, a_start, na_start = optical_fiber_acceptance(radius, acceptance_angle, mid, refracted_2, refracted_endpoints_1, refracted_endpoints_2)
    accepted = np.hstack((a_start, a_end - a_start))  # Combine start and end points for accepted vectors
    not_accepted = np.hstack((na_start, na_end - na_start))  # Combine start and end points for not accepted vectors

    return initial, initial_2, accepted, not_accepted


def point_average_power2(power_absorbed, point, num_vectors, solid_angle, accepted_start, quantum_yield, c, extinction_coeff):
    solid_angle_factor = (len(accepted_start) / num_vectors) * (solid_angle / (4 * np.pi))
    average_vector_power_factor = vectors_average_power2(point, accepted_start, quantum_yield, c, extinction_coeff)
    average_power_factor = power_absorbed * solid_angle_factor * average_vector_power_factor
    return average_power_factor


def point_average_power3(power_absorbed, vector_ratio, solid_angle, quantum_yield):
    average_power_factor = power_absorbed * quantum_yield * vector_ratio * (solid_angle / (4 * np.pi))
    return average_power_factor

def vectors_average_power2(origin, intersection_endpoints, quantum_yield, c, extinction_coeff):
    num_vectors = len(intersection_endpoints)
    if num_vectors > 0:
        vector_powers = np.zeros(num_vectors)
        single_vector_power = 1 / num_vectors * quantum_yield
        for i in range(num_vectors):
            d = distance_3d(origin, intersection_endpoints[i]) / 10
            vector_powers[i] = np.exp(-extinction_coeff * (c * d))
        average_vector_power_factor = np.sum(vector_powers) * single_vector_power
    else:
        average_vector_power_factor = 0
    return average_vector_power_factor

def generate_vectors(h, w, mid, point, num_vectors):
    diff_x = mid[0] - point[0]
    x = point[0] + np.ones(num_vectors) * diff_x
    y = (np.random.rand(num_vectors) * w - w / 2 + mid[1]) - point[1]
    z = (np.random.rand(num_vectors) * h - h / 2 + mid[2]) - point[2]
    vectors = np.column_stack((x, y, z))
    return vectors

def optical_fiber_acceptance(radius, angle, mid, vectors, startpoints, endpoints):
    counter_a = 0
    counter_na = 0
    a_end = np.zeros_like(endpoints)
    na_end = np.zeros_like(endpoints)
    a_start = np.zeros_like(endpoints)
    na_start = np.zeros_like(endpoints)
    for i in range(len(vectors)):
        y = endpoints[i, 1]
        z = endpoints[i, 2]
        y_c = y - mid[1]
        z_c = z - mid[2]
        if np.sqrt(y_c ** 2 + z_c ** 2) > radius:
            na_end[counter_na] = endpoints[i]
            na_start[counter_na] = startpoints[i]
            counter_na += 1
        else:
            if np.arctan(np.abs(vectors[i, 1] / vectors[i, 0])) > angle or np.arctan(np.abs(vectors[i, 2] / vectors[i, 0])) > angle:
                na_end[counter_na] = endpoints[i]
                na_start[counter_na] = startpoints[i]
                counter_na += 1
            else:
                a_end[counter_a] = endpoints[i]
                a_start[counter_a] = startpoints[i]
                counter_a += 1
    a_end = a_end[:counter_a]
    na_end = na_end[:counter_na]
    a_start = a_start[:counter_a]
    na_start = na_start[:counter_na]
    return a_end, na_end, a_start, na_start

def optical_fiber_acceptance_quick(radius, angle, mid, vectors, endpoints):
    a_start = np.zeros(len(endpoints))
    for i in range(len(vectors)):
        y = endpoints[i, 1]
        z = endpoints[i, 2]
        y_c = y - mid[1]
        z_c = z - mid[2]
        if np.sqrt(y_c ** 2 + z_c ** 2) <= radius:
            if np.arctan(np.abs(vectors[i, 1] / vectors[i, 0])) <= angle and np.arctan(np.abs(vectors[i, 2] / vectors[i, 0])) <= angle:
                a_start[i] = 1
    return a_start

def vect_cone_intersection(points, vectors, rho, p1):
    initial_endpoints = np.zeros((len(vectors), 3))
    normal_vectors = np.zeros((len(vectors), 3))
    thet_1 = np.sin(rho)
    thet_2 = np.cos(rho)

    m = 1 / np.tan(rho)
    b = p1[1] - m * p1[0]

    for i in range(len(vectors)):
        m1, b1 = lin_eq([points[i, 0], points[i, 2]], [points[i, 0] + vectors[i, 0], points[i, 2] + vectors[i, 2]])
        r_i = line_intersection(m1, b1, m, b)
        z_f = y_value_at(r_i, m1, b1)

        m2, b2 = lin_eq([points[i, 0], points[i, 1]], [points[i, 0] + vectors[i, 0], points[i, 1] + vectors[i, 1]])
        x_f = intersection_line_with_circle2(m2, b2, r_i)
        y_f = y_value_at(x_f, m2, b2)

        u = (z_f - b) / thet_2
        v = np.arccos(x_f / (u * thet_1))

        initial_endpoints[i, :] = [x_f, y_f, z_f]
        normal_vectors[i, :] = normal_vector_cone(u, v, rho)

    return initial_endpoints, normal_vectors


def normal_vector_cone(u, v, rho):
    a = np.sin(rho)
    b = np.cos(rho)

    return [a * b * u * np.cos(v), a * b * u * np.sin(v), (-(a ** 2)) * u * (np.cos(v)) ** 2 + (-(a ** 2)) * u * (np.sin(v)) ** 2]


def vector_endpoints(vectors, origins, mid_x):
    endpoints = np.zeros((len(vectors), 3))

    for i in range(len(vectors)):
        endpoints[i, :] = vector_endpoint(vectors[i, :], origins[i, :], mid_x)

    return endpoints


def vector_endpoint(vector, origin, mid_x):
    t = (mid_x - origin[0]) / vector[0]
    x_f = mid_x
    y_f = origin[1] + vector[1] * t
    z_f = origin[2] + vector[2] * t

    return [x_f, y_f, z_f]


def refracted_vectors(v_is, norms, n1, n2):
    vectors = np.zeros_like(norms)

    for i in range(len(norms)):
        vectors[i, :] = refracted_vector(v_is[i, :], norms[i, :], n1, n2)

    return vectors


def refracted_vector(v, norm, n1, n2):
    mu = n1 / n2
    dot = v[0] * norm[0] + v[1] * norm[1] + v[2] * norm[2]

    return mu * v + norm * (np.sqrt(1 - (mu ** 2) * (1 - dot ** 2))) - (mu * dot) * norm
