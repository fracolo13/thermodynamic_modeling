# Thermomechanical Simulation 

This repository contains the original MATLAB code of a thermomechanical simulation for slab break-off scenarios.



## Files Description

### MATLAB Files
- `Final_Project_colosimo.m` - **Fixed original code** with matrix singularity solution

## Usage

### MATLAB Version (Fixed)
# Thermomechanical Slab Subduction Simulation (MATLAB)

## Overview
This MATLAB code implements a comprehensive 2D thermomechanical simulation of oceanic slab subduction and potential slab break-off using the marker-in-cell (Eulerian-Lagrangian) method. The simulation couples thermal and mechanical processes to model the evolution of subducting oceanic lithosphere in the mantle.

## Scientific Purpose
The code models key geological processes during subduction:
- **Slab Dynamics**: Evolution of a cold oceanic slab descending into the hot mantle
- **Thermal Evolution**: Heat transfer through conduction, advection, and various heating mechanisms
- **Mechanical Deformation**: Viscous flow governed by temperature and pressure-dependent rheology
- **Slab Break-off**: Potential detachment of the slab due to thermal weakening and gravitational forces

## Physical Domain
- **Grid**: 101×81 Eulerian nodes covering 400×320 km domain
- **Depth Range**: 0-400 km (sticky air layer, lithosphere, and upper mantle)
- **Lagrangian Markers**: 500×400 particles tracking material properties
- **Time Scale**: Geological time steps (~10¹¹-10¹² seconds) over multiple million years

## Key Physical Components

### 1. Material Layers
- **Sticky Air (0-50 km)**: Low-density layer representing atmosphere/surface
  - Density: 1 kg/m³, Viscosity: 10¹⁸ Pa·s, Temperature: 273 K
- **Oceanic Plate/Slab (50-100 km + subducted portion)**: Cold, dense oceanic lithosphere
  - Density: 3400 kg/m³, Viscosity: 10²³ Pa·s, Variable temperature with depth
- **Mantle (>100 km)**: Hot, convecting mantle material  
  - Density: 3250 kg/m³, Viscosity: 10²⁰ Pa·s, Temperature: 1573 K
- **Necking Zone**: Weakened region in slab with reduced yield strength for potential break-off

### 2. Slab Geometry
The code initializes a tilted slab geometry:
- **Surface Plate**: Horizontal oceanic lithosphere at the surface
- **Subducting Slab**: Inclined slab extending from 100-250 km depth
- **Temperature Gradient**: Linear temperature increase with depth in the slab
- **Necking Area**: Localized weakening zone for modeling slab detachment

### 3. Governing Equations

#### Mechanical (Stokes Flow)
```
∇·σ + ρg = 0    (Momentum conservation)
∇·v = 0         (Mass conservation - incompressible flow)
```
Where σ is the stress tensor, ρ is density, g is gravity, and v is velocity.

#### Thermal Evolution  
```
ρCₚ ∂T/∂t = ∇·(k∇T) + H_rad + H_s + H_a
```
Where:
- **H_rad**: Radioactive heating (2×10⁻⁸ W/m³ in plates, 3×10⁻⁸ W/m³ in mantle)
- **H_s**: Shear heating from viscous dissipation
- **H_a**: Adiabatic heating from pressure changes

### 4. Advanced Numerical Methods

#### Non-linear Rheology
- **Temperature-dependent viscosity**: η = η₀ exp(E/(RT)) for different rock types
- **Pressure-dependent**: Includes pressure effects on activation energy
- **Yield stress**: Plastic yielding when stress exceeds material strength
- **Strain rate dependence**: Power-law relationship for creep deformation

#### Marker-in-Cell Method
- **Property Interpolation**: Weighted averaging from markers to Eulerian grid
- **Advection**: 4th-order Runge-Kutta scheme for accurate particle tracking
- **Heat Conservation**: Conservative temperature interpolation maintaining energy balance

#### Thermal Conductivity
Temperature-dependent thermal conductivity:
```
k(T) = k₀ + a/(T + b)
```
Where k₀, a, b are material-specific constants.

### 5. Boundary Conditions

#### Mechanical
- **Top/Bottom**: Free slip (tangential stress = 0)
- **Sides**: No normal flow, free tangential slip
- **Internal**: Continuity of velocity and stress

#### Thermal  
- **Top**: Fixed temperature (273 K - surface conditions)
- **Bottom**: Fixed temperature (1573 K - mantle conditions)  
- **Sides**: Symmetry conditions (zero heat flux)

### 6. Heating Mechanisms

#### Adiabatic Heating
```
H_a = T·α·ρ·g·v_y
```
Heating due to pressure changes during vertical motion.

#### Shear Heating
```
H_s = σ_xx·ε_xx + σ_yy·ε_yy + 2σ_xy·ε_xy
```
Viscous dissipation from deformation.

#### Subgrid Diffusion
Advanced thermal diffusion accounting for temperature variations smaller than grid resolution.

### 7. Adaptive Time Stepping
The code employs sophisticated time step control:
- **CFL Condition**: dt < min(dx/v_max, dy/v_max) for numerical stability
- **Thermal Stability**: Ensures proper thermal evolution
- **Mechanical Convergence**: Iterative solution until velocity convergence

## Simulation Workflow

1. **Initialization**: Set up geometry, material properties, and initial conditions
2. **Time Loop**: For each geological time step:
   - Interpolate marker properties to Eulerian grid
   - Assemble and solve mechanical equations (Stokes flow)
   - Calculate thermal source terms (adiabatic, shear heating)
   - Solve thermal diffusion equation
   - Update material properties based on new T,P conditions
   - Advect markers using 4th-order Runge-Kutta
   - Apply subgrid thermal diffusion
   - Advance to next time step

## Expected Results
The simulation produces:
- **Velocity Fields**: Flow patterns showing slab descent and mantle circulation
- **Temperature Evolution**: Thermal structure showing slab cooling and mantle heating  
- **Stress/Strain Fields**: Deformation patterns and stress concentrations
- **Potential Slab Break-off**: Detachment when thermal weakening reduces slab strength

## Applications
This model is used to study:
- Subduction zone thermal structure and dynamics
- Slab stagnation and break-off processes  
- Mantle flow patterns around subducting slabs
- Heat transfer in convergent plate boundaries
- Temporal evolution of subduction systems

## Technical Notes
- Uses MATLAB's sparse matrix solvers for efficient computation
- Employs marker-in-cell method for accurate material tracking
- Implements advanced thermal diffusion with subgrid corrections
- Includes comprehensive rheological models for realistic deformation
- Matrix regularization prevents numerical issues with poorly conditioned systems

The code represents a state-of-the-art implementation of coupled thermomechanical modeling for understanding fundamental processes in subduction zone dynamics and slab evolution.

---

## References

This code implements numerical methods for:
- Finite difference discretization of Stokes equations
- Thermal diffusion with advection
- Marker-in-cell particle methods
- Sparse linear algebra and iterative solvers

The physical model is suitable for studying geodynamic processes like slab subduction and break-off events.
