# Lunar Surface Mission - Return Phases
**Role:** Team Leader (4-person engineering team)  
**Domain:** Flight Dynamics, GNC, Astrodynamics, Space Robotics  

## Overview
This repository contains the MATLAB framework developed for designing and simulating the orbital return phases of a Lunar Module (LM) to a Mothership (MS). The project spans from two-body Hohmann transfers and high-fidelity perturbed orbital dynamics (J2, 3rd-body) to Circular Restricted Three-Body Problem (CR3BP) targeting and 5-DOF robotic berthing operations.

## Core Modules & Engineering Methodologies

### 1. Orbital Transfer & Phasing (`main_part1.m`)
* **Analytical Modeling:** Computed delta-V budgets and phasing requirements for a LM-MS Hohmann transfer.
* **Numerical Validation:** Validated analytical Keplerian solutions (Newton-Raphson) against high-precision numerical integration (`ode45` with 1e-14 tolerances).

### 2. Relative Dynamics & Proximity Operations (`main_part2.m`)
* **Linear vs. Nonlinear Dynamics:** Modeled the free-drift and forced V-bar approach in the LVLH frame using the Hill-Clohessy-Wiltshire (HCW) State Transition Matrix.
* **Targeting Problem:** Solved a non-linear two-impulse targeting maneuver using `fsolve` to bring the LM to a 50m hold point, validating the limits of the linear HCW approximation.

### 3. Perturbed Orbital Dynamics (`main_part3.m`)
* **Cowell's Method:** Integrated the $J_2$ effect and Earth's third-body gravitational perturbations directly into the equations of motion.
* **Osculating Elements:** Extracted and visualized the secular drift in semi-major axis, eccentricity, and orbital angles over time.

### 4. CR3BP Transfer & Shooting Method (`main_bonus_CR3BP.m`)
* **High-Fidelity Targeting:** Implemented a numerical Shooting Method (`fsolve` + `ode45`) within the non-dimensionalized Earth-Moon CR3BP to compensate for third-body radial and phase errors in the standard Hohmann transfer.
* **Numerical Stability:** Monitored the Jacobi Constant ($C$) throughout the integration to guarantee energy conservation within machine precision bounds.
![CR3BP Transfer](./images/Hohmann_Transfer.png)

### 5. Robotic Berthing Dynamics (`main_part4.m`)
* **Spatial Kinematics:** Modeled a 5R Mothership robotic arm using Homogeneous Transformation Matrices (HTM) and analytically derived the Geometric and System Jacobians.
* **Inverse Kinematics & Trajectory:** Solved IK near singularities using the Damped Least Squares (DLS) method and planned smooth joint-space trajectories via cubic polynomial interpolation.

## Prerequisites & Installation
To run the simulations, the following MATLAB Toolboxes are required:
* Optimization Toolbox (for `fsolve` targeting)
* Symbolic Math Toolbox (for analytical Jacobians)
* Robotics System Toolbox

## Usage
Clone the repository and run any of the main module scripts (e.g., `main_part2.m`). Each script automatically clears the workspace, runs the specific mission phase simulation, and plots the relevant astrodynamics data.
