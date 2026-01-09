from typing import Callable
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
        self._e: float = 0.6 # eccentricity
        self._rp: float = 2.0  # semimajor axis
        self._update_orbital_params_erp(self._e, self._rp)
        self._update_orbit_type(self._e)


    @property
    def e(self) -> float:
        return self._e

    @e.setter
    def e(self, e: str) -> None:
        e = float(e)
        self._e = e
        self._update_orbital_params_erp(e, self._rp)
        self._update_orbit_type(e)

    @property
    def rp(self) -> float:
        return self._rp

    @rp.setter
    def rp(self, rp: str) -> None:
        rp = float(rp)
        self._rp = rp
        self._update_orbital_params_erp(self._e, rp)

    @property
    def a(self) -> float:
        return self._a

    @a.setter
    def a(self, a: str) -> None:
        a = float(a)
        self._a = float(a)

    @property
    def b(self) -> float:
        return self._b

    @b.setter
    def b(self, b: float | int) -> None:
        self._b = float(b)

    @property
    def p(self) -> float:
        return self._p

    @p.setter
    def p(self, value: float) -> None:
        self._p = float(value)


    @property
    def ra(self) -> float:
        return self._ra

    @ra.setter
    def ra(self, value: float | int) -> None:
        self._ra = float(value)

    @property
    def orbit_eq(self) -> PerifocalOrbitEq:
        return PerifocalOrbitEq(self._e, self._p)

    def orbital_angles(self):
        if self._e <= 1:
            return np.linspace(0, 2*pi, 1000)

        delta = 0.0001
        return np.linspace(-self._asymptote_anomaly + delta, self._asymptote_anomaly - delta, 1000)

    def _update_orbital_params_erp(self, e: float, rp: float):
        self._p: float = self._orbital_param_erp(e, rp)
        self._a: float = self._semimajor_axis_erp(e, rp)
        self._b: float = self._semiminor_axis_erp(e, rp)
        self._ra: float = self._apoapsis_ep(e, self._p)
        self._asymptote_anomaly: float = self._asymptote_anomaly_e(e)
        print(f"p = {self._p}")
        print(f"a = {self._a}")
        print(f"b = {self._b}")
        print(f"ra = {self._ra}")

    def _update_orbit_type(self, e: float) -> None:
        orbit_types: dict[str, bool] = {
            "_circular" : e == 0,
            "_elliptical": 0 < e < 1,
            "_parabolic": e == 1,
            "_hyperbolic": e > 1
        }
        for type, val in list(orbit_types.items()):
            self.__setattr__(type, val)

    def _orbital_param_erp(self, e: float, rp: float) -> float:
        return rp*(1 + e)

    def _semimajor_axis_erp(self, e: float, rp: float) -> float:
        return rp/(1 - e)

    def _semiminor_axis_erp(self, e: float, rp: float) -> float:
        if e > 1:
            return rp*np.sqrt(e**2 - 1)/(1 - e)

        return rp*np.sqrt(1 - e**2)/(1 - e)

    def _apoapsis_ep(self, e: float, p: float) -> float:
        return p*(1/(1 - e))

    def _asymptote_anomaly_e(self, e: float) -> float:
        """Calculate the true anomaly of the asymptote for hyperbolic orbits using the eccentricity"""
        if e > 1:
            return np.arccos(-1/e)

        return np.nan

