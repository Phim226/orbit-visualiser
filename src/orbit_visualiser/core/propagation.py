from scipy.integrate import solve_ivp
from math import ceil
import numpy as np
from numpy.typing import NDArray
from orbit_visualiser.core.orbit import Orbit

# TODO: Implement CLI input.

def two_body_ode(t: float, state: NDArray[np.float64], mu: float) -> NDArray[np.float64]:
    """
    The second order ordinary differential equation describing the relative motion of a body in a
    two body problem.

    Parameters
    ----------
    t : float
        Unused time argument. Necessary for the way scipy uses functions to numerically integrate.
    state : NDArray[np.float64]
        The current state [x, dx], with x in km and dx in km/s.
    mu : float
        The gravitational parameter of the system. For satellites around a planet we use the
        gravitational parameter of the planet.

    Returns
    -------
    NDArray[np.float64]
        The values to be numerically integrated [dx, d^2x]
    """
    pos = state[:3]
    r = np.linalg.norm(pos)

    vel = state[3:]
    acc: NDArray[np.float64] = -(mu/r**3)*pos

    return np.concatenate((vel, acc))

def get_init_conditions_from_orbit(orbit: Orbit) -> NDArray[np.float64]:
    """
    Takes orbit object and returns the initial conditions (position and velocity) for orbit propagation.

    Parameters
    ----------
    orbit: Orbit
        Orbit object to get the initial conditions from

    Returns
    -------
    NDArray[np.float64]
        The initial conditions [r, v] in the form of a concatenated numpy array.
    """
    return np.concatenate((orbit.position, orbit.velocity))

def run_orbit_prop(
        orbit: Orbit, t_end: float, adaptive_step_size: bool = True, period_frac_per_step: int = 500
) -> list[tuple[float, NDArray[np.float64]]]:
    """
    Run the orbit propagation for the satellite. Uses the RK45 algorithm.

    Parameters
    ----------
    orbit : Orbit
        The orbit on which the satellite is being propagated.
    t_end : float
        The end time of the propagation. It should be on the order of at most 10 orbital periods
        since the integrator (RK45) isn't symplectic, so will suffer from energy drift over long
        propagations.
    adaptive_step_size: bool, optional
        Boolean flag deciding whether to use adaptive step sizes. If set to false then the step size
        is determined by the period_frac_per_step optional argument, by default True.
    period_frac_per_step : int, optional
        The time step size as a fraction of the orbital period. So with the default value of 500, each
        time step is 1/500 of the orbital period, by default 500.

    Returns
    -------
    list[tuple[float, NDArray[np.float64]]]
        List with each element containing propagation time and the propagated state vector at that
        time. Extracting the last element of this list will give the state vector at the end of
        the propagation interval.
    """
    if not adaptive_step_size:
        step_size = max(100, ceil((t_end/orbit.orbital_period)*period_frac_per_step))
        t = np.linspace(0, t_end, step_size)
    else:
        t = None

    init_conditions = get_init_conditions_from_orbit(orbit)

    result = solve_ivp(
        two_body_ode,
        [0, t_end],
        init_conditions,
        args = (orbit.mu,),
        t_eval = t,
        rtol=1e-10,
        atol=1e-12
    )

    print(result.message)
    return [(t, y) for t, y in zip(result.t, result.y.T)]

if __name__ == "__main__":

    orbit = Orbit.from_orbital_elements(e = 0.0, rp = 50_000.0, nu = 0.0, raan = 0.0, i = 0.0,
                                        omega = 0.0, mu = 398_600.0,)

    sol = run_orbit_prop(orbit, orbit.orbital_period)

    init_conditions = get_init_conditions_from_orbit(orbit)
    print(init_conditions)
    print(sol.y[:, -1])

    r0 = np.linalg.norm(init_conditions[:3])
    rf = np.linalg.norm(sol.y[:3, -1])

    print(f"Difference in radius after propagating: {rf - r0}")
