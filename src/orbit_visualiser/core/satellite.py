import numpy as np
from orbit_visualiser.core import Orbit, CentralBody

# TODO: Write methods for flight angle, perifocal positions and characteristic energy.
# TODO: Refactor so that nu is in radians instead of degrees.
class Satellite():


    def __init__(self, orbit: Orbit, central_body: CentralBody, nu: float = 0.00):
        self._orbit = orbit
        self._central_body = central_body
        self._nu: float = nu # true anomaly in °

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

    def update_satellite_properties(self) -> None:
        nu, mu, e, rp = self._nu, self._central_body.mu, self._orbit.e, self._orbit.rp
        h, t_asymp = self._specific_ang_momentum(mu, rp, e), self._orbit.t_asymp
        self._h = h
        self._x, self._y = self._position(nu, t_asymp)
        self._v_azim = self._azimuthal_velocity(mu, h, e, nu, t_asymp)
        self._v_radial = self._radial_velocity(mu, h, e, nu)
        self._v = self._velocity(self._v_azim, self._v_radial)
        self._r = self._radius(h, mu, e, nu, t_asymp)
        self._eps = self._specific_energy(mu, h, e)
        self._v_esc = self._escape_velocity(nu, mu, self._r, t_asymp)
        self._v_inf = self._excess_velocity(mu, abs(self._orbit.a))

    @staticmethod
    def _specific_ang_momentum(mu: float, rp: float, e: float) -> float:
        return np.sqrt(mu*rp*(1 + e))

    def _position(self, nu: float, t_asymp: float) -> tuple[float]:
        if np.isclose(abs(nu), t_asymp, atol = 0.005, rtol = 0):
                return np.inf, np.inf

        orbit_eq = self._orbit.orbit_eq
        return orbit_eq.x(np.deg2rad(nu)), orbit_eq.y(np.deg2rad(nu))


    def _azimuthal_velocity(self, mu: float, h: float, e: float, nu: float, t_asymp: float) -> float:
        if np.isclose(abs(nu), t_asymp, atol = 0.005, rtol = 0):
            return 0.0

        return (mu/h)*(1 + e*np.cos(np.deg2rad(nu)))

    def _radial_velocity(self, mu: float, h: float, e: float, nu: float) -> float:
        if self._orbit.orbit_type == "circular":
            return 0.0

        return (mu/h)*e*np.sin(np.deg2rad(nu))

    def _velocity(self, v_azim: float, v_radial: float) -> float:
        return np.sqrt(v_azim**2 + v_radial**2)

    def _radius(self, h: float, mu: float, e: float, nu: float, t_asymp: float) -> float:
        if np.isclose(abs(nu), t_asymp, atol = 0.005, rtol = 0):
            return np.inf

        return (h**2/mu)/(1 + e*np.cos(np.deg2rad(nu)))

    def _specific_energy(self, mu: float, h: float, e: float) -> float:
        if self._orbit.orbit_type == "parabolic":
            return 0.0

        return -0.5*(mu/h)**2*(1 - e**2)


    def _escape_velocity(self, nu: float, mu: float, r: float, t_asymp: float) -> float:
        if np.isclose(abs(nu), t_asymp, atol = 0.005, rtol = 0):
            return 0.0

        return np.sqrt(2*mu/r)

    def _excess_velocity(self, mu: float, a: float) -> float:
        if self._orbit.is_closed:
            return np.nan

        return np.sqrt(mu/a)
