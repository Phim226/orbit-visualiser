# Orbit Visualiser

The Orbit Visualiser is a 2D Keplerian orbit visualisation tool for modelling the motion of a satellite around a central body, written in Python. The focus of this tool is on orbital geometry and satellite kinematics. Has basic orbital propagation.

## Features

Orbits are modelled and visualised in the perifocal frame (the central body remains fixed at the origin of the coordinate frame), with the orbital geometry parametrised using:
  - Eccentricity
  - Radius of periapsis

The central body has a radius of 6738km and an adjustable gravitational parameter. The radius of periapsis has therefore been given a lower limit of 6739km (meaning that we are assuming this Earth sized body has no atmosphere, or at least that any atmosphere has no effect on orbital motion). The true anomaly of the orbiting satellite is also adjustable to evaluate the kinematic state at different orbital positions. Various orbital and kinematic quantities are
calculated and displayed, including:
  - Semi-major/minor axis
  - Radius of apoapsis
  - True anomaly of the asymptote (for parabolic and hyperbolic trajectories)
  - Specific angular momentum
  - Radial velocity
  - Azimuthal velocity
  - Escape velocity
  - Hyperbolic excess velocity
  - Time since periapsis

## Notes on the physical model

The orbiting satellite is assumed to have negligible mass. All higher-order perturbations, such as atmospheric drag or oblateness, are ignored.

## Potential future improvements

- Introduce argument of periapsis as an orbital geometry parameter
- Improve simulation of orbital motion
- Expand orbital modelling to 3 dimensions
- Include perturbations into the model
- Expand modelling to R3BP or N-body system
- Visualise gravitational field

## Requirements

- Python 3.8 or later
- Git (optional but recommended, for cloning the repository)

## Python dependencies (also listed in requirements.txt)

- numpy
- matplotlib
- scipy

## How to use
Clone the repository:

```bash
git clone https://github.com/Phim226/orbit-visualiser.git
```

Navigate to the project folder:

```bash
cd orbit-visualiser
```
Create a virtual environment:

```python
python3 -m venv .venv
```

Activate the environment:

### Windows
*Windows Powershell*:

```powershell
.venv\Scripts\Activate.ps1
```
- **Warning**: On some Windows systems, Powershell might block this command depending on the script execution policy. If you see an error like:

  ```powershell
  ...\.venv\Scripts\Activate.ps1 cannot be loaded because running scripts is disabled on this system.
  ```
  Then you can run:

  ```powershell
  Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope Process
  ```
  to **temporarily** allow scripts in the current terminal session.

*Windows Command Prompt*:

```bash
.venv\Scripts\activate.bat
```
This should work without changing any execution policies.
### Linux/macOS

```bash
source .venv/bin/activate
```

If you see (.venv) at the start of your prompt, then the virtual environment is active.

Install dependencies and the package in editable mode (make sure this is done within the virtual environment):

```python
python3 -m pip install --upgrade pip
pip install -r requirements.txt
```

The ```-e``` in requirements.txt ensures that the package is editable, which ensures that any changes to the source code are immediately reflected when running the program.

Run the program:

```python
python src/orbit_visualiser/main.py
```

When you're finished using the program or making changes then run:

```bash
deactivate
```
to deactivate the virtual environment. When returning to the program, activate venv again using the commands above. So long as the project folder or .venv folder hasn't changed then you shouldn't need to reinstall any dependencies.
