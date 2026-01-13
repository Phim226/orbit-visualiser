# Orbit Visualiser

The Orbit Visualiser is a 2D Keplerian orbit visualisation tool for modelling the orbit of a satellite around a central body, written in Python. The focus of this tool is on orbital geometry and satellite kinematics. Currently there is no simulation-based functionality.

## Features

Orbits are modelled and visualised in the perifocal frame (the central body remains fixed at the origin of the coordinate frame), with the orbital geometry parametrised using:
  - Eccentricity
  - Radius of periapsis

The central body has an adjustable gravitational parameter. The true anomaly of the orbiting satellite is also adjustable to evaluate the kinematic state at different orbital positions. Various orbital properties and kinematics quantities are 
calculated and displayed, including:
  - Semi-major/minor axis
  - Radius of apoapsis
  - True anomaly of the asymptote (for parabolic and hyperbolic trajectories)
  - Specific angular momentum
  - Radial velocity
  - Azimuthal velocity
  - Escape velocity
  - Hyperbolic excess velocity

## Notes on the physical model

The orbiting satellite is assumed to have negligible mass. All higher order perturbations, such as atmospheric drag or oblateness, are ignored.

## Potential future improvements

- Calculate time since periapsis
- Simulate orbital motion
- Include perturbations into the model
- Expand modelling to R3BP or N-body system
- Visualise gravitational field

## Python dependencies

- Numpy
- Matplotlib
