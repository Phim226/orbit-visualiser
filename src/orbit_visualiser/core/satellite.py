from orbit_visualiser.core import Orbit

# TODO: Write methods for calculating satellite state (velocity, distance from pericenter etc)
class Satellite():

    def __init__(self, orbit: Orbit):
        self._orbit = orbit
        self._nu: float = 0.0 # true anomaly in °

    @property
    def nu(self) -> float:
        return self._nu

    @nu.setter
    def nu(self, nu: float) -> None:
        self._nu = nu