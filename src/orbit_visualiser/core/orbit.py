from dataclasses import dataclass
from math import pi
import numpy as np

# TODO: Split formulae from Orbit class.
# TODO: Update orbit properties when eccentricity of radius of periapsis is set rather than calling it explicitly outside the class
class Orbit():


    def __init__(self, e: float = 0.0, rp: float = 50_000):
        self._e: float = e # Eccentricity
        self._rp: float = rp # Radius of periapsis in km
        self.update_orbital_properties()


    @property
    def e(self) -> float:
        return self._e

    @e.setter
    def e(self, e: float) -> None:
        self._e = e

    @property
    def rp(self) -> float:
        return self._rp

    @rp.setter
    def rp(self, rp: float) -> None:
        self._rp = rp

    @property
    def a(self) -> float:
        return self._a

    @property
    def b(self) -> float:
        return self._b

    @property
    def p(self) -> float:
        return self._p

    @property
    def ra(self) -> float:
        return self._ra

    @property
    def t_asymp(self) -> float:
        return self._t_asymp

    @property
    def turn_angle(self) -> float:
        return self._turn_angle

    @property
    def aim_rad(self) -> float:
        return self._aim_rad

    @property
    def orbit_type(self) -> str:
        return self._orbit_type

    @property
    def is_closed(self) -> bool:
        return self._is_closed

    def orbital_angles(self) -> tuple[float]:
        if self._e < 1:
            return 0, 2*pi

        return -self._t_asymp, self._t_asymp

    def update_orbital_properties(self):
        e, rp = self._e, self._rp

        self._update_orbit_type(e)

        self._p: float = self._orbital_param_erp(e, rp)
        a = self._semimajor_axis_erp(e, rp)
        self._a: float = a
        b = self._semiminor_axis_erp(e, a)
        self._b: float = b
        self._ra: float = self._apoapsis_erp(e, a)
        self._t_asymp: float = self._asymptote_anomaly_e(e)
        self._turn_angle: float = self._turning_angle_e(e)
        self._aim_rad: float = self._aiming_radius_erp(b)

    def _update_orbit_type(self, e: float) -> None:
        """Determine orbit type and if it's closed using the eccentricity"""
        if e == 0:
            self._orbit_type = "circular"
            self._is_closed = True

        elif 0 < e < 1:
            self._orbit_type = "elliptical"
            self._is_closed = True

        elif e == 1:
            self._orbit_type = "parabolic"
            self._is_closed = False

        else:
            self._orbit_type = "hyperbolic"
            self._is_closed = False

    def _orbital_param_erp(self, e: float, rp: float) -> float:
        """Calculate the orbital parameter p using the eccentricity and radius of periapsis"""
        return rp*(1 + e)

    def _semimajor_axis_erp(self, e: float, rp: float) -> float:
        """Calculate the semimajor axis a using the eccentricity and radius of periapsis"""
        if self._orbit_type != "parabolic":
            return rp/(1 - e)

        return np.inf

    def _semiminor_axis_erp(self, e: float, a: float) -> float:
        """Calculate the semiminor axis b using the eccentricity and semi-major axis"""
        if self._orbit_type == "hyperbolic":
            return a*np.sqrt(e**2 - 1)

        elif self._is_closed:
            return a*np.sqrt(1 - e**2)

        return np.inf

    def _apoapsis_erp(self, e: float, a: float) -> float:
        """Calculate the radius of apoapsis ra using the eccentricity and semi-major axis"""
        if self._orbit_type != "parabolic":
            return a*(1 + e)

        return np.inf

    def _asymptote_anomaly_e(self, e: float) -> float:
        """Calculate the true anomaly of the asymptote for open orbits using the eccentricity"""
        if not self._is_closed:
            return np.arccos(-1/e)

        return np.nan

    def _turning_angle_e(self, e: float) -> float:
        """Calculate the turning angle for open orbits using the eccentricity"""
        if not self._is_closed:
            return 2*np.arcsin(1/e)

        return np.nan

    def _aiming_radius_erp(self, b: float) -> float:
        """Calculate the aiming radius for hyperbolic orbits using the semi-minor axis"""
        if self._orbit_type == "hyperbolic":
            return -b

        return np.nan

@dataclass
class CentralBody():
        mu: float = 398600 # gravitational parameter in km³/s² = Gm
        r: float = 6378 # radius in km
