from scipy.integrate import solve_ivp
from functools import partial
import numpy as np
from numpy.typing import NDArray
from orbit_visualiser.core.orbit import Orbit
from orbit_visualiser.core.satellite import Satellite, CentralBody


def two_body_pf_ode(mu: float, t: float, state: NDArray[np.float64], ) -> NDArray[np.float64]:
    x, y, v_x, v_y = state

    r = np.hypot(x, y)

    a_x = -(mu/r**3)*x
    a_y = -(mu/r**3)*y

    return np.array([v_x, v_y, a_x, a_y])

def get_init_conditions(satellite: Satellite, at: str = "periapsis") -> NDArray[np.float64]:
    if at == "periapsis":
        sat_copy = Satellite(satellite._orbit, satellite._central_body)

        return np.concatenate((sat_copy.pos_pf, sat_copy.vel_pf))

    return np.concatenate((satellite.pos_pf, satellite.vel_pf))

def run_orbit_prop(satellite: Satellite, init_conditions: NDArray[np.float64], t_span: tuple[float], period_frac_per_step: int = 500):
    t = np.linspace(0, satellite.period, period_frac_per_step)

    sol = solve_ivp(
        partial(two_body_pf_ode, satellite._central_body.mu),
        t_span,
        init_conditions,
        t_eval = t,
        rtol=1e-10,
        atol=1e-12
    )
    return sol

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
