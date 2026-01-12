from orbit_visualiser.core import Orbit

class Satellite():

    def __init__(self, orbit: Orbit):
        self._orbit = Orbit
        self._nu: float = 0.0 # true anomaly in °

    @property
    def nu(self) -> float:
        return self._nu

    @nu.setter
    def nu(self, nu: float) -> None:
        self._nu = nu