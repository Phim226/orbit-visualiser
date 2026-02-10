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

sat = Satellite(Orbit(), CentralBody())

r = sat.r
mu = sat._central_body.mu
period = sat.period

x = r
y = 0
v_x = 0
v_y = np.sqrt(mu/r)

t = np.linspace(0, period, 1000)

sol = solve_ivp(partial(two_body_pf_ode, mu), [0, sat.period], [x, y, v_x, v_y], t_eval = t, rtol=1e-10, atol=1e-12)

print(sol.y[:, -1])