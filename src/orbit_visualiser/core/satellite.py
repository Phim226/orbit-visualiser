from orbit_visualiser.core import Orbit, CentralBody

# TODO: Write methods for calculating satellite state (velocity, distance from pericenter etc)
class Satellite():

    def __init__(self, orbit: Orbit, central_body: CentralBody):
        self._orbit = orbit
        self._central_body = central_body
        self._nu: float = 0.0 # true anomaly in °

        self._update_satellite_properties()

    @property
    def nu(self) -> float:
        return self._nu

    @nu.setter
    def nu(self, nu: float) -> None:
        self._nu = nu

    def _update_satellite_properties(self) -> None:
        pass