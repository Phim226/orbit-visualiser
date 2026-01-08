from typing import Callable
import numpy as np


class PerifocalOrbitEq():

    def __init__(self, e: float, a: float):
        p: float = a*(1 - e**2) # orbital parameter
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
        self._e = float(value)

    @property
    def a(self) -> float:
        return self._a

    @a.setter
    def a(self, value: float | int) -> None:
        self._a = float(value)

    @property
    def b(self) -> float:
        return self._b

    @b.setter
    def b(self, value: float | int) -> None:
        self._b = float(value)

    def orbit_eq(self) -> PerifocalOrbitEq:
        return PerifocalOrbitEq(self._e, self._a)
