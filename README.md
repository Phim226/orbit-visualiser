# Orbit Visualiser

The Orbit Visualiser is a 2D Keplerian orbit visualisation tool for modelling the motion of a satellite around a central body, written in Python. The focus of this tool is on orbital geometry and satellite kinematics. Currently there is no simulation-based functionality.

## Features

Orbits are modelled and visualised in the perifocal frame (the central body remains fixed at the origin of the coordinate frame), with the orbital geometry parametrised using:
  - Eccentricity
  - Radius of periapsis

The central body has an adjustable gravitational parameter. The true anomaly of the orbiting satellite is also adjustable to evaluate the kinematic state at different orbital positions. Various orbital and kinematic quantities are 
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

The orbiting satellite is assumed to have negligible mass. All higher-order perturbations, such as atmospheric drag or oblateness, are ignored.

## Potential future improvements

- Introduce argument of periapsis as an orbital geometry parameter
- Calculate time since periapsis
- Simulate orbital motion
- Include perturbations into the model
- Expand modelling to R3BP or N-body system
- Visualise gravitational field

## Requirements

- Python 3.8 or later
- Git (optional but recommended, for cloning the repository)

## Python dependencies 

- numpy
- matplotlib

These are also listed in the requirements.txt file.

## How to use
Clone the repository:

```bash
git clone https://github.com/Phim226/orbit-visualiser.git
```

Navigate to project folder:

```bash
cd orbit-visualiser
```

Install dependencies:

```bash
pip install -r requirements.txt
```

Run the program:

```bash
python src/orbit_visualiser/main.py
```
