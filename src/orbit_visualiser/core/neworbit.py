from dataclasses import dataclass
from typing import Sequence
import numpy as np
from numpy.typing import NDArray
from orbit_visualiser.core.astrodynamics.keplerian.state import state_pf_from_e_rp

@dataclass
class CentralBody():
    """
    Represents the body around which a satellite orbits. The default values represent Earth.

    Parameters
    ----------
    mu : float
        The gravitational parameter (km^3/s^2), default = 398600.0
    r : float
        The radius (km), default = 6378.0
    """
    mu: float = 398600.0
    r: float = 6378.0


class NewOrbit():
    """
    Represents the instantaneous analytical orbit of a satellite. If there are no perturbations
    (oblateness, drag, a third body, etc) then for given inputs the orbit will be fixed. If
    perturbations are present then that causes orbital elements to change over time, represented
    by the changing state of this object.

    Parameters
    ----------
    r : Sequence | NDArray[np.float64]
        The perifocal position vector of the satellite (km)
    v : Sequence | NDArray[np.float64]
        The perifocal velocity vector of the satellite (km/s)
    mu : float
        The gravitational parameter of the central body (km^3/s^2)

    Other Constructors
    ----------
    Class can be constructed from orbital elements. Use Orbit.from_orbital_elements.

    Parameters
    ----------
    e : float
        The eccentricity
    rp : float
        The radius of periapsis (km)
    mu : float
        The gravitational parameter of the central body (km^3/s^2)
    nu : float
        The true anomaly of the satellite (rads)
    """

    def __init__(
            self,
            r: Sequence | NDArray[np.float64],
            v: Sequence | NDArray[np.float64],
            mu: float
    ):
        self._r = r
        self._v = v
        self._mu = mu


    @classmethod
    def from_orbital_elements(cls, e: float, rp: float, mu: float, nu: float):
        r, v = state_pf_from_e_rp(e, rp, mu, nu)
        return cls(r, v, mu)

    def compute_orbital_elements(
            self,
            r: Sequence | NDArray[np.float64],
            v: Sequence | NDArray[np.float64],
            mu: float
    ) -> None:
        self.e # use eccentricity vector equation

