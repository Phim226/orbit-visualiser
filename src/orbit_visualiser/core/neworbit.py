from dataclasses import dataclass
from functools import cached_property
from typing import Sequence
import numpy as np
from numpy.typing import NDArray
from orbit_visualiser.core.astrodynamics.keplerian.state import state_pf_from_e_rp
from orbit_visualiser.core.astrodynamics.keplerian.elements import (eccentricity_from_state, semi_parameter_from_momentum,
                                                                    periapsis, semimajor_axis, semiminor_axis,
                                                                    apoapsis, asymptote_anomaly, turning_angle,
                                                                    aiming_radius, orbital_period, mean_motion)
from orbit_visualiser.core.astrodynamics.keplerian.dynamics import (specific_orbital_energy, characteristic_energy,
                                                                    excess_velocity, specific_ang_momentum_from_state)
from orbit_visualiser.core.astrodynamics.keplerian.classification import orbit_type
from orbit_visualiser.core.astrodynamics.types import OrbitType

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

@dataclass
class NewOrbit():
    """
    Represents the instantaneous analytical orbit of a satellite. If there are no perturbations
    (oblateness, drag, a third body, etc) then for given inputs the orbit will be fixed. If
    perturbations are present then that causes orbital elements to change over time, represented
    by the changing state of this object.

    Parameters
    ----------
    position : Sequence | NDArray[np.float64]
        The perifocal position vector of the satellite (km)
    velocity : Sequence | NDArray[np.float64]
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
    position : Sequence | NDArray[np.float64]
    velocity : Sequence | NDArray[np.float64]
    mu : float

    @cached_property
    def eccentricity(self) -> float:
        return eccentricity_from_state(self.position, self.velocity, self.mu)

    @cached_property
    def orbit_type(self) -> OrbitType:
        return orbit_type(self.eccentricity)

    @cached_property
    def semi_parameter(self) -> OrbitType:
        return semi_parameter_from_momentum(specific_ang_momentum_from_state(self.position, self.velocity), self.mu)

    @cached_property
    def radius_of_periapsis(self) -> OrbitType:
        return periapsis(self.semi_parameter, self.eccentricity)

    @cached_property
    def semimajor_axis(self) -> OrbitType:
        return semimajor_axis(self.orbit_type, self.eccentricity, self.radius_of_periapsis)

    @cached_property
    def semiminor_axis(self) -> OrbitType:
        return semiminor_axis(self.orbit_type, self.eccentricity, self.semimajor_axis)

    @cached_property
    def radius_of_apoapsis(self) -> OrbitType:
        return apoapsis(self.orbit_type, self.eccentricity, self.semimajor_axis)

    @cached_property
    def asymptote_anomaly(self) -> OrbitType:
        return asymptote_anomaly(self.orbit_type, self.eccentricity)

    @cached_property
    def turning_angle(self) -> OrbitType:
        return turning_angle(self.orbit_type, self.eccentricity)

    @cached_property
    def aiming_radius(self) -> OrbitType:
        return aiming_radius(self.orbit_type, self.semiminor_axis)

    @cached_property
    def orbital_period(self) -> OrbitType:
        return orbital_period(self.orbit_type, self.mu, self.semimajor_axis)

    @cached_property
    def mean_motion(self) -> OrbitType:
        return mean_motion(self.orbit_type, self.orbital_period, self.mu, self.semi_parameter, self.semimajor_axis)

    @cached_property
    def specific_energy(self) -> OrbitType:
        return specific_orbital_energy(self.orbit_type, self.mu, self.semimajor_axis)

    @cached_property
    def characteristic_energy(self) -> OrbitType:
        return characteristic_energy(self.orbit_type, self.mu, self.semimajor_axis)

    @cached_property
    def hyperbolic_excess_velocity(self) -> OrbitType:
        return excess_velocity(self.orbit_type, self.mu, self.semimajor_axis)

    @classmethod
    def from_orbital_elements(cls, e: float, rp: float, mu: float, nu: float, asymptote_anomaly: float):
        """
        Alternative constructor for the Orbit class. Takes the orbital elements eccentricity,
        radius of perapsis, the gravitational parameter and the true anomaly as arguments.

        Parameters
        ----------
        e : float
            Eccentricity
        rp : float
            Radius of periapsis (km)
        mu : float
            The gravitational parameter of the central body (km^3/s^2)
        nu : float
            The true anomaly of the satellite (rads)

        Returns
        -------
        Orbit
            A new instance of Orbit
        """
        r, v = state_pf_from_e_rp(e, rp, mu, nu, asymptote_anomaly)
        return cls(r, v, mu)


