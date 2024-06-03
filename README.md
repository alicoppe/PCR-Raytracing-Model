# PCR-Raytracing-Model

## Introduction

This simulation creates a parametric model for a PCR tube containing a variable amount of liquid. It aims to assess the proportion of fluorescent light emitted by fluorophores in the solution that reaches an optical fiber for detection. The PCR tube is modeled as a two-part system consisting of an exterior shell and an interior shell containing the fluid. Measurements are made in millimeters and are based on actual data.

## Model Description

### Parametric Model

The PCR tube is modeled as a cylinder with concentric washers representing the liquid volume. Each washer has a specified height \(dz\) along the z-axis. The emission of light from each fluorophore is modeled using the Beer-Lambert law, and the rays are traced using a 3D version of Snell's Law.

### Mathematical Formulation

#### Parametric Equations

1. **Cylinder (PCR Tube):**
   - Radius: $\( r \)$
   - Height: $\( h \)$
   - Equation: 
     $$x = r \cos(\theta), \quad y = r \sin(\theta), \quad z = z$$

2. **Washer (Slice of Liquid):**
   - Inner Radius: $\( r_1 \)$
   - Outer Radius: $\( r_2 \)$
   - Height: $\( dz \)$
   - Volume of washer:
     $$V_{\text{washer}} = \pi (r_2^2 - r_1^2) \cdot dz$$

3. **Solid Angle (\(\Omega\)) in Spherical Coordinates:**
   - $$\Omega = \int_0^{2\pi} \int_0^{\theta} \sin(\phi) \, d\phi \, d\theta$$

#### Beer-Lambert Law

The Beer-Lambert Law describes the absorption of light as it travels through a medium:
$$I(z) = I_0 e^{-\alpha z)$$
Where:
- $\( I(z) \)$ is the intensity of light at depth $\( z \)$
- $\( I_0 \)$ is the initial intensity of light
- $\( \alpha \)$ is the absorption coefficient

#### Solid Angle and Spherical Coordinates

To calculate the solid angle subtended by the optical fiber from a point in the liquid:
$$\Omega = \int_0^{2\pi} \int_0^{\theta} \sin(\phi) \, d\phi \, d\theta$$

#### Snell's Law in 3D

Snell's Law describes the refraction of light at the interface between two media:
$$n_1 \sin(\theta_1) = n_2 \sin(\theta_2)$$
Where:
- $\( n_1 \)$ and $\( n_2 \)$ are the refractive indices of the two media
- $\( \theta_1 \)$ is the angle of incidence
- $\( \theta_2 \)$ is the angle of refraction

In 3D, Snell's Law can be applied to each component of the incident ray.


## Simulation Steps

### Volume Discretization

1. Divide the PCR tube into concentric washers along the z-axis.
2. Each washer is discretized into a uniform distribution of points.

### Ray Tracing

1. For each point in a washer, generate rays within the solid angle subtended by the optical fiber.
2. Apply Snell's Law to determine the direction and refraction of each ray.
3. Calculate the proportion of rays that reach the optical fiber.

### Power Calculation

1. Use the Beer-Lambert Law to calculate the fluorescence intensity at each point.
2. Integrate the contributions of all points within a washer to find the average power reaching the optical fiber.
3. Sum the power contributions from all washers to get the total power.

## Dependencies

- Python 3.10.10
- NumPy
- Matplotlib (for visualization)
