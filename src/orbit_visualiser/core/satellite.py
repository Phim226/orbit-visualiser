import numpy as np
from orbit_visualiser.core import Orbit, CentralBody

# TODO: Write methods for calculating satellite state (velocity, distance from pericenter etc)
class Satellite():

    def __init__(self, orbit: Orbit, central_body: CentralBody):
        self._orbit = orbit
        self._central_body = central_body
        self._nu: float = 0.0 # true anomaly in °

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
    def r(self) -> float:
        return self._r

    @property
    def eps(self) -> float:
        return self._eps

    def update_satellite_properties(self) -> None:
        nu, mu, e, rp = np.deg2rad(self._nu), self._central_body.mu, self._orbit.e, self._orbit.rp
        h = self._specific_ang_momentum(mu, rp, e)
        self._h = h
        self._v_azim = self._azimuthal_velocity(mu, h, e, nu)
        self._v_radial = self._radial_velocity(mu, h, e, nu)
        self._v = self._velocity(self._v_azim, self._v_radial)
        self._r = self._radius(h, mu, e, nu)
        self._eps = self._specific_energy(mu, h, e)

    @staticmethod
    def _specific_ang_momentum(mu: float, rp: float, e: float) -> float:
        return np.sqrt(mu*rp*(1 + e))

    @staticmethod
    def _azimuthal_velocity(mu: float, h: float, e: float, nu: float) -> float:
        return (mu/h)*(1 + e*np.cos(nu))

    @staticmethod
    def _radial_velocity(mu: float, h: float, e: float, nu: float) -> float:
        return (mu/h)*e*np.sin(nu)

    @staticmethod
    def _velocity(v_azim: float, v_radial: float) -> float:
        return np.sqrt(v_azim**2 + v_radial**2)

    @staticmethod
    def _radius(h: float, mu: float, e: float, nu: float) -> float:
        if np.isclose(nu, -np.pi) or np.isclose(nu, np.pi):
            return np.inf
        return (h**2/mu)/(1 + e*np.cos(nu))

    @staticmethod
    def _specific_energy(mu: float, h: float, e: float) -> float:
        return -0.5*(mu/h)**2*(1 - e**2)
