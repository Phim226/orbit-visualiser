from scipy.integrate import solve_ivp
from functools import partial
from math import ceil
import numpy as np
from numpy.typing import NDArray
from orbit_visualiser.core.orbit import Orbit
from orbit_visualiser.core.satellite import Satellite, CentralBody


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
        The current state [x, dx]

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

def get_init_conditions(satellite: Satellite, at: str = "periapsis") -> NDArray[np.float64]:
    """
    Takes a Satellite object and returns the initial conditions (position and velocity) for orbit propagation.

    Parameters
    ----------
    satellite : Satellite
        The satellite object, which contains information about its analytical orbit.
    at : str, optional
        An optional parameter to choose where the in the analytical orbit the initial conditions are
        taken from, by default "periapsis".

    Returns
    -------
    NDArray[np.float64]
        The initial conditions in the form of a concatenated numpy array.
    """
    if at == "periapsis":
        sat_copy = Satellite(satellite._orbit, satellite._central_body)

        return np.concatenate((sat_copy.pos_pf, sat_copy.vel_pf))

    return np.concatenate((satellite.pos_pf, satellite.vel_pf))

def run_orbit_prop(satellite: Satellite, init_conditions: NDArray[np.float64], t_end: float, period_frac_per_step: int = 500):
    """
    Run the orbit propagation for the satellite. Uses the RK45 algorithm.

    Parameters
    ----------
    satellite : Satellite
        The satellite object to be propagated.
    init_conditions : NDArray[np.float64]
        The initial conditions for the propagation.
    t_end : float
        The end time of the propagation. It should be on the order of a single orbital period
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
    t = np.linspace(0, t_end, ceil((t_end/satellite.period)*period_frac_per_step))

    result = solve_ivp(
        partial(two_body_pf_ode, satellite._central_body.mu),
        [0, t_end],
        init_conditions,
        t_eval = t,
        rtol=1e-10,
        atol=1e-12
    )
    return result

if __name__ == "__main__":

    satellite = Satellite(Orbit(), CentralBody())
    init_conditions = get_init_conditions(satellite)
    t_span = [0, satellite.period]

    sol = run_orbit_prop(satellite, init_conditions, t_span)
    print(init_conditions)
    print(sol.y[:, -1])

    r0 = np.linalg.norm(init_conditions[:2])
    rf = np.linalg.norm(sol.y[:2, -1])

    print(f"Difference in radius after propagating: {rf - r0}")
