from typing import Callable
from dataclasses import dataclass
from math import pi
import numpy as np


class PerifocalOrbitEq():

    def __init__(self, e: float, p: float):
        self._x: Callable[[float], float] = lambda t : p*(np.cos(t)/(1+e*np.cos(t)))
        self._y: Callable[[float], float] = lambda t : p*(np.sin(t)/(1+e*np.cos(t)))

    @property
    def x(self) -> Callable[[float], float]:
        return self._x

    @property
    def y(self) -> Callable[[float], float]:
        return self._y

class Orbit():

    def __init__(self):
        self._e: float = 0.6
        self._rp: float = 10_000
        self._update_orbital_params_erp(self._e, self._rp)
        self._update_orbit_type(self._e)


    @property
    def e(self) -> float:
        return self._e

    @e.setter
    def e(self, e: float) -> None:
        self._e = e
        self._update_orbital_params_erp(e, self._rp)
        self._update_orbit_type(e)

    @property
    def rp(self) -> float:
        return self._rp

    @rp.setter
    def rp(self, rp: float) -> None:
        self._rp = rp
        self._update_orbital_params_erp(self._e, rp)

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
    def t_asymp(self, units: str = "degrees") -> float:
        if units == "degrees":
            return np.degrees(self._t_asymp)

        return self._t_asymp

    @property
    def turn_angle(self, units: str = "degrees") -> float:
        if units == "degrees":
            return np.degrees(self._turn_angle)

        return self._turn_angle

    @property
    def aim_rad(self) -> float:
        return self._aim_rad

    @property
    def orbit_eq(self) -> PerifocalOrbitEq:
        return PerifocalOrbitEq(self._e, self._p)

    def orbital_angles(self):
        if self._e < 1:
            return np.linspace(0, 2*pi, 1000)

        delta = 0.0001
        return np.linspace(-self._t_asymp + delta, self._t_asymp - delta, 1000)

    def _update_orbital_params_erp(self, e: float, rp: float):
        self._p: float = self._orbital_param_erp(e, rp)
        self._a: float = self._semimajor_axis_erp(e, rp)
        self._b: float = self._semiminor_axis_erp(e, rp)
        self._ra: float = self._apoapsis_erp(e, rp)
        self._t_asymp: float = self._asymptote_anomaly_e(e)
        self._turn_angle: float = self._turning_angle_e(e)
        self._aim_rad: float = self._aiming_radius_erp(e, rp)

    def _update_orbit_type(self, e: float) -> None:
        if e == 0:
            self._orbit_type = "circular"

        elif 0 < e < 1:
            self._orbit_type = "elliptical"

        elif e == 1:
            self._orbit_type = "parabolic"

        else:
            self._orbit_type = "hyperbolic"

    def _orbital_param_erp(self, e: float, rp: float) -> float:
        """Calculate the orbital parameter p using the eccentricity and radius of periapsis"""
        return rp*(1 + e)

    def _semimajor_axis_erp(self, e: float, rp: float) -> float:
        """Calculate the semimajor axis a using the eccentricity and radius of periapsis"""
        if e != 1:
            return rp/(1 - e)

        return np.inf

    def _semiminor_axis_erp(self, e: float, rp: float) -> float:
        """Calculate the semiminor axis b using the eccentricity and radius of periapsis"""
        if e > 1:
            return rp*np.sqrt(e**2 - 1)/(1 - e)

        elif e < 1:
            return rp*np.sqrt(1 - e**2)/(1 - e)

        return np.inf

    def _apoapsis_erp(self, e: float, rp: float) -> float:
        """Calculate the radius of apoapsis ra using the eccentricity and radius of periapsis"""
        if e != 1:
            return rp*((1 + e)/(1 - e))

        return np.inf

    def _asymptote_anomaly_e(self, e: float) -> float:
        """Calculate the true anomaly of the asymptote for hyperbolic orbits using the eccentricity"""
        if e >= 1:
            return np.arccos(-1/e)

        return np.nan

    def _turning_angle_e(self, e: float) -> float:
        """Calculate the turning angle for hyperbolic orbits using the eccentricity"""
        if e >= 1:
            return 2*np.arcsin(1/e)

        return np.nan

    def _aiming_radius_erp(self, e: float, rp: float) -> float:
        if e > 1:
            return rp*np.sqrt((e + 1)/(e - 1))

        return np.nan

@dataclass
class CentralBody():
        mu = 398600 # gravitational parameter in km³/s² = Gm
        r = 6378 # radius in km
