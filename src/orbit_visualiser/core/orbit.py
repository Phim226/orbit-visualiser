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
        self._a: float = 5.0  # semimajor axis
        self._update_orbital_params_ea(self._e, self._a)


    @property
    def e(self) -> float:
        return self._e

    @e.setter
    def e(self, e: str) -> None:
        e = float(e)
        if e > 1:
            self._asymptote_anomaly = np.arccos(-1/e)
            if self._a > 0 :
                self._a = -np.abs(self._a)
                self._b = -np.abs(self._b)

        if (e < 1 and self._a < 0):
            self._a = np.abs(self._a)
            self._b = np.abs(-self._b)

        self._e = e
        self._update_orbital_params_ea(e, self._a)

    @property
    def a(self) -> float:
        return self._a

    @a.setter
    def a(self, a: str) -> None:
        a = float(a)
        self._a = a
        self._update_orbital_params_ea(self._e, a)

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
    def rp(self) -> float:
        return self._rp

    @rp.setter
    def rp(self, value: float | int) -> None:
        self._rp = float(value)

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

    def _update_orbital_params_ea(self, e: float, a: float):
        self._b: float = self._semiminor_axis_ea(e, a) # semiminor axis
        self._p: float = self._orbital_param_ea(e, a) # orbital parameter
        self._rp: float = self._periapsis_ep(e, self._p) # radius of periapsis
        self._ra: float = self._apoapsis_ep(e, self._p) # radius of apoapsis
        print(f"b = {self._b}")
        print(f"p = {self._p}")
        print(f"rp = {self._rp}")
        print(f"ra = {self._ra}")

    def _semiminor_axis_ea(self, e: float, a: float) -> float:
        if e > 1:
            return a*np.sqrt(e**2 - 1)

        return a*np.sqrt(1 - e**2)

    def _orbital_param_ea(self, e: float, a: float) -> float:
        return a*(1 - e**2)

    def _periapsis_ep(self, e: float, p: float) -> float:
        return p*(1/(1 + e))

    def _apoapsis_ep(self, e: float, p: float) -> float:
        return p*(1/(1 - e))
