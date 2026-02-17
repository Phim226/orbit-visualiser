from scipy.integrate import solve_ivp
from functools import partial
from math import ceil
import numpy as np
from numpy.typing import NDArray
from orbit_visualiser.core.orbit import Orbit

# TODO: Implement CLI input.

def two_body_pf_ode(mu: float, t: float, state: NDArray[np.float64], ) -> NDArray[np.float64]:
    """
    The second order ordinary differential equation describing the relative motion of a body in a
    two body problem.

    Parameters
    ----------
    mu : float
        The gravitational parameter of the system. For satellites around a planet we use the
        gravitational parameter of the planet.
    t : float
        Unused time argument. Necessary for the way scipy uses functions to numerically integrate.
    state : NDArray[np.float64]
        The current state [x, dx], with x in km and dx in km/s.

    Returns
    -------
    NDArray[np.float64]
        The values to be numerically integrated [dx, d^2x]
    """
    x, y, v_x, v_y = state

    r = np.hypot(x, y)

    a_x = -(mu/r**3)*x
    a_y = -(mu/r**3)*y

    return np.array([v_x, v_y, a_x, a_y])

def get_init_conditions_from_orbit(orbit: Orbit) -> NDArray[np.float64]:
    """
    Takes orbital elements and returns the initial conditions (position and velocity) for orbit propagation.

    Parameters
    ----------
    orbit: NewOrbit
        Orbit object to get the initial conditions from


    Returns
    -------
    NDArray[np.float64]
        The initial conditions [r, v] in the form of a concatenated numpy array.
    """
    return np.concatenate((orbit.position, orbit.velocity))

def run_orbit_prop(orbit: Orbit, t_end: float, period_frac_per_step: int = 500):
    """
    Run the orbit propagation for the satellite. Uses the RK45 algorithm.

    Parameters
    ----------
    orbit: NewOrbit
        The orbit on which the satellite is being propagated.
    t_end : float
        The end time of the propagation. It should be on the order of at most 10 orbital periods
        since the integrator (RK45) isn't symplectic, so will suffer from energy drift over long
        propagations.
    period_frac_per_step : int, optional
        The time step size as a fraction of the orbital period. So with the default value of 500, each
        time step is 1/500 of the orbital period, by default 500.

    Returns
    -------
    OdeResult (see scipy documentation)
        The object containing the results of the numerical integration. The array containing the
        [x, dx] values can be accessed by calling result.y
    """
    t = np.linspace(0, t_end, ceil((t_end/orbit.orbital_period)*period_frac_per_step))

    init_conditions = get_init_conditions_from_orbit(orbit)

    result = solve_ivp(
        partial(two_body_pf_ode, orbit.mu),
        [0, t_end],
        init_conditions,
        t_eval = t,
        rtol=1e-10,
        atol=1e-12
    )
    return result

if __name__ == "__main__":

    orbit = Orbit.from_orbital_elements(e = 0.0, rp = 50_000.0, mu = 398_600.0, nu = 0.0)

    sol = run_orbit_prop(orbit, orbit.orbital_period)

    init_conditions = get_init_conditions_from_orbit(orbit)
    print(init_conditions)
    print(sol.y[:, -1])

    r0 = np.linalg.norm(init_conditions[:2])
    rf = np.linalg.norm(sol.y[:2, -1])

    print(f"Difference in radius after propagating: {rf - r0}")
