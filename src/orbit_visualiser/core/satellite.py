import numpy as np
from math import pi
from orbit_visualiser.core import Orbit, CentralBody

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

    def update_satellite_properties(self) -> None:
        nu, mu, e, rp = self._nu, self._central_body.mu, self._orbit.e, self._orbit.rp
        h, t_asymp, a = self._specific_ang_momentum(mu, rp, e), self._orbit.t_asymp, abs(self._orbit.a)
        self._h = h
        self._x, self._y = self._position(nu, t_asymp)
        self._v_azim = self._azimuthal_velocity(mu, h, e, nu, t_asymp)
        self._v_radial = self._radial_velocity(mu, h, e, nu)
        self._v = self._velocity(mu, h, e, nu)
        self._r = self._radius(h, mu, e, nu, t_asymp)
        self._eps = self._specific_energy(mu, h, e)
        self._v_esc = self._escape_velocity(nu, mu, self._r, t_asymp)
        self._v_inf = self._excess_velocity(mu, a)
        self._gam = self._flight_angle(e, nu, t_asymp)
        self._c3 = self._characteristic_energy(mu, a)
        self._t = self._orbital_period(mu, a)

    @staticmethod
    def _specific_ang_momentum(mu: float, rp: float, e: float) -> float:
        return np.sqrt(mu*rp*(1 + e))

    def _position(self, nu: float, t_asymp: float) -> tuple[float]:
        if np.isclose(abs(nu), t_asymp, atol = 0.0001, rtol = 0):
                return np.inf, np.inf

        orbit_eq = self._orbit.orbit_eq
        return orbit_eq.x(nu), orbit_eq.y(nu)


    def _azimuthal_velocity(self, mu: float, h: float, e: float, nu: float, t_asymp: float) -> float:
        if np.isclose(abs(nu), t_asymp, atol = 0.0001, rtol = 0):
            return 0.0

        return (mu/h)*(1 + e*np.cos(nu))

    def _radial_velocity(self, mu: float, h: float, e: float, nu: float) -> float:
        if self._orbit.orbit_type == "circular":
            return 0.0

        return (mu/h)*e*np.sin(nu)

    def _velocity(self, mu: float, h: float, e: float, nu: float) -> float:
        return (mu/h)*np.sqrt(e**2 + 2*e*np.cos(nu) + 1)

    def _radius(self, h: float, mu: float, e: float, nu: float, t_asymp: float) -> float:
        if np.isclose(abs(nu), t_asymp, atol = 0.0001, rtol = 0):
            return np.inf

        return (h**2/mu)/(1 + e*np.cos(nu))

    def _specific_energy(self, mu: float, h: float, e: float) -> float:
        if self._orbit.orbit_type == "parabolic":
            return 0.0

        return -0.5*(mu/h)**2*(1 - e**2)

    def _escape_velocity(self, nu: float, mu: float, r: float, t_asymp: float) -> float:
        if np.isclose(abs(nu), t_asymp, atol = 0.0001, rtol = 0):
            return 0.0

        return np.sqrt(2*mu/r)

    def _excess_velocity(self, mu: float, a: float) -> float:
        if self._orbit.is_closed:
            return np.nan

        return np.sqrt(mu/a)

    def _flight_angle(self, e: float, nu: float, t_asymp: float) -> float:
        if np.isclose(nu, t_asymp, atol = 0.0001, rtol = 0):
            return pi/2

        elif np.isclose(nu, -t_asymp, atol = 0.0001, rtol = 0):
            return -pi/2

        return np.arctan((e*np.sin(nu))/(1 + e*np.cos(nu)))

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
