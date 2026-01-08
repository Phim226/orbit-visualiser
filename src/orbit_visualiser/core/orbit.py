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
        self._rp: float = 2.0 # radius of periapsis
        self._ra: float = 8.0 # radius of apoapsis
        self._e: float = 0.6 # eccentricity
        self._a: float = 5.0  # semimajor axis
        self._b: float = 4  # semiminor axis

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
    def e(self) -> float:
        return self._e

    @e.setter
    def e(self, value: float | int) -> None:
        e = float(value)
        if e > 1:
            self._asymptote_anomaly = np.arccos(-1/e)
            if self._a > 0 :
                self._a = -self._a
                self._b = -self._b

        if (e < 1 and self._a < 0):
            self._a = -self._a
            self._b = -self._b

        self._e = e
        self.p = self._orbital_param_ea(e, self._a)

    @property
    def a(self) -> float:
        return self._a

    @a.setter
    def a(self, a: float) -> None:
        a = float(a)
        self._a = a
        self.p = self._orbital_param_ea(self._e, a)

    @property
    def b(self) -> float:
        return self._b

    @b.setter
    def b(self, value: float | int) -> None:
        self._b = float(value)

    @property
    def p(self) -> float:
        return self._p

    @p.setter
    def p(self, value: float) -> None:
        self._p = float(value)

    @property
    def orbit_eq(self) -> PerifocalOrbitEq:
        return PerifocalOrbitEq(self._e, self._p)

    def orbital_angles(self):
        if self._e <= 1:
            return np.linspace(0, 2*pi, 1000)

        delta = 0.0001
        return np.linspace(-self._asymptote_anomaly + delta, self._asymptote_anomaly - delta, 1000)
