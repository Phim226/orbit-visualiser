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

## Planned Features

Some GUI elements such as additional display options are currently commented out because they are under development. These will be implemented in future updates.

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

## Python dependencies

- numpy
- matplotlib
- scipy

These are also listed in requirements.txt.

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

**Note**: The current GUI layout hasn't been optimised for Linux, so there is overlapping of certain GUI elements due to differences in widget rendering and layout handling.

## Using Orbit Propagation
As of the current version (v0.4.5) there is no CLI or GUI functionality for the orbit propagation. Navigate to src/orbit_visualiser/core/propagation.py, and go to the bottom of the file to the code snippet:
```python
if __name__ == "__main__":

    orbit = Orbit.from_orbital_elements(e = 0.0, rp = 50_000.0, mu = 398_600.0, nu = 0.0)

    sol = run_orbit_prop(orbit, orbit.orbital_period)

    init_conditions = get_init_conditions_from_orbit(orbit)
    print(init_conditions)
    print(sol.y[:, -1])

    r0 = np.linalg.norm(init_conditions[:2])
    rf = np.linalg.norm(sol.y[:2, -1])

    print(f"Difference in radius after propagating: {rf - r0}")
```
The orbit propagation function takes the Orbit object and propagation end time (the start time is always 0, and begins the propagation at the true anomaly (nu) set in the constructor Orbit.from_orbital_elements). The values set in the code snippet above result in the propagation of a circular orbit or radius 50,000km, starting at the perifocal position (50000, 0) lasting for one orbital period. For non-circular orbits nu = 0.0 is periapsis. The true anomaly is in radians, so nu = pi is apoapsis.

In order to get accurate results you should keep the end time on the order of 10 or at most 100 periods. The integrators used by scipy aren't symplectic, so there is significant energy drift over longer propagations.

There is a third optional argument of run_orbit_prop called period_frac_per_step: int = 500, which determines the time step as a function of the fraction of a single period. So the default value results in a time step equal to 1/500 of the orbital period (rounded up to the nearest integer).

The current print statements are just quick sanity checks to show that after a whole number of orbits the satellite is returning to roughly the same spot and with roughly the same velocity.

As an example if you wanted to propagate the ISS then you might edit the above code snippet to:
```python
if __name__ == "__main__":
    # Approximate eccentricity and radius of periapsis (in km) for the ISS, with mu representing
    # the gravitational parameter of earth in km^3/s^2
    iss_orbit = Orbit.from_orbital_elements(e = 0.0002267, rp = 6778, mu = 398_600.0, nu = 0.0)

    # Runs the propagation for 1 period, returning the ISS back to periapsis
    sol = run_orbit_prop(iss, iss.orbital_period)

    init_conditions = get_init_conditions_from_orbit(iss_orbit)
    print(init_conditions)

    # Grabs the array of solutions (sol.y) and prints the last entry
    print(sol.y[:, -1])

    r0 = np.linalg.norm(init_conditions[:2])
    rf = np.linalg.norm(sol.y[:2, -1])

    print(f"Difference in radius after propagating: {rf - r0}")
```
The running the script will run the propagation and print out the results.

Alternatively, you can pass the initial position, velocity and gravitational parameter directly to Orbit. By using
```python
orbit = Orbit(position = np.array([50_000.0, 0.0]), velocity = np.array([0.0, 2.82]), mu = 398600.0)
```
Since we are in the perifocal frame then all orbits are oriented with the periapsis and apoapsis on the x-axis, so the true anomaly is calculated from
```python
np.atan2(0.0, 50_000.0) = 0.0
```
You can pass this Orbit object into run_orbit_prop (as well as your desired run time) just as before.