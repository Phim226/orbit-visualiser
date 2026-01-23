import numpy as np
from math import pi
from orbit_visualiser.core import Orbit, CentralBody

# TODO: Write equations for eccentric anomaly of hyperbolic trajectories.
# TODO: Write equations for mean anomaly of open trajectories.
# TODO: Write equations for time since periapsis for elliptical and open trajectories.
class Satellite():


    def __init__(self, orbit: Orbit, central_body: CentralBody, nu: float = 0.00):
        self._orbit = orbit
        self._central_body = central_body
        self._nu: float = nu # true anomaly in rads

        self.update_satellite_properties()

    @property
    def nu(self) -> float:
        return self._nu

    @nu.setter
    def nu(self, nu: float) -> None:
        self._nu = nu

    @property
    def h(self) -> float:
        return self._h

    @property
    def v_azim(self) -> float:
        return self._v_azim

    @property
    def v_radial(self) -> float:
        return self._v_radial

    @property
    def v(self) -> float:
        return self._v

    @property
    def v_esc(self) -> float:
        return self._v_esc

    @property
    def v_inf(self) -> float:
        return self._v_inf

    @property
    def r(self) -> float:
        return self._r

    @property
    def eps(self) -> float:
        return self._eps

    @property
    def x(self) -> float:
        return self._x

    @property
    def y(self) -> float:
        return self._y

    @property
    def gam(self) -> float:
        return self._gam

    @property
    def c3(self) -> float:
        return self._c3

    @property
    def t(self) -> float:
        return self._t

    @property
    def n(self) ->float:
        return self._n

    @property
    def t_p(self) -> float:
        return self._t_p

    def update_satellite_properties(self) -> None:
        nu, mu, e, p = self._nu, self._central_body.mu, self._orbit.e, self._orbit.p
        self._h, t_asymp, a = self._specific_ang_momentum(mu, p), self._orbit.t_asymp, self._orbit.a

        # Helper quantities
        den = 1 + e*np.cos(nu)
        mu_over_h = mu/self._h

        # Geometry
        self._x, self._y = self._position(nu, t_asymp)
        self._r = self._radius(nu, t_asymp, p, den)

        # Kinematics
        v_azim = self._azimuthal_velocity(nu, t_asymp, mu_over_h, den)
        self._v_azim = v_azim
        v_radial = self._radial_velocity(mu_over_h, e, nu)
        self._v_radial = v_radial
        self._v = self._velocity(v_azim, v_radial)
        self._v_esc = self._escape_velocity(nu, mu, self._r, t_asymp)
        self._v_inf = self._excess_velocity(mu, abs(a))
        self._gam = self._flight_angle(e, nu, t_asymp, den)

        # Energy
        self._eps = self._specific_energy(mu, a)
        self._c3 = self._characteristic_energy(mu, abs(a))

        # Time
        self._t = self._orbital_period(mu, a)
        self._n = self._mean_motion(self._t)
        self._t_p = self._time_since_periapsis(nu, self._t)

    @staticmethod
    def _specific_ang_momentum(mu: float, p: float) -> float:
        return np.sqrt(mu*p)

    def _position(self, nu: float, t_asymp: float) -> tuple[float]:
        if np.isclose(abs(nu), t_asymp, atol = 0.0001, rtol = 0):
                return np.inf, np.inf

        orbit_eq = self._orbit.orbit_eq
        return orbit_eq.x(nu), orbit_eq.y(nu)


    def _azimuthal_velocity(self, nu: float, t_asymp: float, mu_over_h: float, den: float) -> float:
        if np.isclose(abs(nu), t_asymp, atol = 0.0001, rtol = 0):
            return 0.0

        return mu_over_h*den

    def _radial_velocity(self, mu_over_h: float, e: float, nu: float) -> float:
        if self._orbit.orbit_type == "circular":
            return 0.0

        return mu_over_h*e*np.sin(nu)

    def _velocity(self, v_azim: float, v_radial: float) -> float:
        return np.hypot(v_azim, v_radial)

    def _radius(self, nu: float, t_asymp: float, p: float, den: float) -> float:
        if np.isclose(abs(nu), t_asymp, atol = 0.0001, rtol = 0):
            return np.inf

        return p/den

    def _specific_energy(self, mu: float, a: float) -> float:
        if self._orbit.orbit_type == "parabolic":
            return 0.0

        return -mu/(2*a)

    def _escape_velocity(self, nu: float, mu: float, r: float, t_asymp: float) -> float:
        if np.isclose(abs(nu), t_asymp, atol = 0.0001, rtol = 0):
            return 0.0

        return np.sqrt(2*mu/r)

    def _excess_velocity(self, mu: float, a: float) -> float:
        if self._orbit.is_closed:
            return np.nan

        return np.sqrt(mu/a)

    def _flight_angle(self, e: float, nu: float, t_asymp: float, den: float) -> float:
        if np.isclose(nu, t_asymp, atol = 0.0001, rtol = 0):
            return pi/2

        elif np.isclose(nu, -t_asymp, atol = 0.0001, rtol = 0):
            return -pi/2

        return np.arctan2(e*np.sin(nu), den)

    def _characteristic_energy(self, mu: float, a: float) -> float:
        if self._orbit.orbit_type == "parabolic":
            return 0.0

        c3 = mu/a
        if self._orbit.is_closed:
            c3 *= -1

        return c3

    def _orbital_period(self, mu: float, a: float) -> float:
        if self._orbit.is_closed:
            return (2*pi/np.sqrt(mu))*np.sqrt(a)**3

        return np.nan

    def _mean_motion(self, t: float) -> float:
        if self._orbit.is_closed:
            return 2*pi/t

        return np.nan

    def _eccentric_anomaly(self, e: float, nu: float):
        orbit_type = self._orbit.orbit_type
        if orbit_type == "circular":
            return nu

        elif orbit_type == "elliptical":
            return 2*np.arctan(np.sqrt((1- e)/(1 + e))*np.tan(nu/2))

        elif orbit_type == "hyperbolic":
            return np.nan

        return np.nan

    def _mean_anomaly(self, e: float, nu: float, E: float):
        orbit_type = self._orbit.orbit_type
        if orbit_type == "circular":
            return nu

        elif orbit_type == "elliptical":
            return E - e*np.sin(E)

        elif orbit_type == "parabolic":
            return np.nan

        elif orbit_type == "hyperbolic":
            return np.nan

    def _time_since_periapsis(self, nu: float, t: float) -> float:
        if self._orbit.orbit_type == "circular":
            return t*nu/(2*pi)

        return np.nan