from dataclasses import dataclass
from functools import cached_property
from typing import Sequence
import numpy as np
from numpy.typing import NDArray
from orbit_visualiser.core.astrodynamics.keplerian.state import state_pf_from_e_rp
from orbit_visualiser.core.astrodynamics.keplerian.elements import (eccentricity_from_state, semi_parameter_from_momentum,
                                                                    radius_of_periapsis, semimajor_axis, semiminor_axis,
                                                                    radius_of_apoapsis, asymptote_anomaly, turning_angle,
                                                                    aiming_radius, orbital_period, mean_motion, inclination,
                                                                    node_line, right_ascen_of_ascending_node, argument_of_periapsis,
                                                                    eccentricity_vector_from_state)
from orbit_visualiser.core.astrodynamics.keplerian.dynamics import (specific_orbital_energy, characteristic_energy,
                                                                    excess_speed, specific_ang_momentum_from_state)
from orbit_visualiser.core.astrodynamics.keplerian.classification import orbit_type, orbit_motion_type
from orbit_visualiser.core.astrodynamics.transformations import perifocal_to_eci_trans_mat
from orbit_visualiser.core.astrodynamics.types import OrbitType, OrbitMotion

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
class Orbit():
    """
    Represents the instantaneous analytical orbit of a satellite. If there are no perturbations
    (oblateness, drag, a third body, etc) then for given inputs the orbit will be fixed. If
    perturbations are present then that causes orbital elements to change over time, represented
    by the changing state of this object.

    Parameters
    ----------
    position : Sequence | NDArray[np.float64]
        The ECI position vector of the satellite (km)
    velocity : Sequence | NDArray[np.float64]
        The ECI velocity vector of the satellite (km/s)
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
        The true anomaly of the satellite (rad)
    raan : float
        The right ascension of the ascending node (rad), default = 0.0
    i : float
        The orbital inclination (rad), default = 0.0
    omega : float
        The argument of periapsis (rad), default = 0.0
    """
    position : Sequence | NDArray[np.float64]
    velocity : Sequence | NDArray[np.float64]
    mu : float

    def __post_init__(self):
        if not np.all(np.isfinite(self.position)) or not np.all(np.isfinite(self.velocity)):
            raise ValueError("State vectors must only contain finite values.")

    # ------ Core orbital elements -------
    @cached_property
    def specific_angular_momentum(self) -> NDArray[np.float64]:
        return specific_ang_momentum_from_state(self.position, self.velocity)

    @cached_property
    def eccentricity(self) -> float:
        return eccentricity_from_state(self.position, self.velocity, self.mu)

    @cached_property
    def inclination(self) -> float:
        return inclination(self.specific_angular_momentum)

    @cached_property
    def right_ascen_of_ascend_node(self) -> float:
        return right_ascen_of_ascending_node(node_line(self.specific_angular_momentum))

    @cached_property
    def argument_of_periapsis(self) -> float:
        return argument_of_periapsis(node_line(self.specific_angular_momentum), self.eccentricity)

    # ------------------------------------

    @cached_property
    def eccentricity_vector(self) -> NDArray[np.float64]:
        return eccentricity_vector_from_state(self.position, self.velocity, self.mu)

    @cached_property
    def orbit_type(self) -> OrbitType:
        return orbit_type(self.eccentricity)

    @cached_property
    def orbit_motion_type(self) -> OrbitMotion:
        return orbit_motion_type(self.inclination)

    @cached_property
    def semi_parameter(self) -> float:
        return semi_parameter_from_momentum(np.linalg.norm(self.specific_angular_momentum), self.mu)

    @cached_property
    def radius_of_periapsis(self) -> float:
        return radius_of_periapsis(self.semi_parameter, self.eccentricity)

    @cached_property
    def semimajor_axis(self) -> float:
        return semimajor_axis(self.eccentricity, self.radius_of_periapsis)

    @cached_property
    def semiminor_axis(self) -> float:
        return semiminor_axis(self.eccentricity, self.semimajor_axis)

    @cached_property
    def radius_of_apoapsis(self) -> float:
        return radius_of_apoapsis(self.eccentricity, self.semimajor_axis)

    @cached_property
    def asymptote_anomaly(self) -> float:
        return asymptote_anomaly(self.eccentricity)

    @cached_property
    def turning_angle(self) -> float:
        return turning_angle(self.eccentricity)

    @cached_property
    def aiming_radius(self) -> float:
        return aiming_radius(self.semiminor_axis)

    @cached_property
    def orbital_period(self) -> float:
        return orbital_period(self.eccentricity, self.mu, self.semimajor_axis)

    @cached_property
    def mean_motion(self) -> float:
        return mean_motion(self.eccentricity, self.mu, self.semi_parameter, self.semimajor_axis)

    @cached_property
    def specific_energy(self) -> float:
        return specific_orbital_energy(self.mu, self.semimajor_axis)

    @cached_property
    def characteristic_energy(self) -> float:
        return characteristic_energy(self.mu, self.semimajor_axis)

    @cached_property
    def hyperbolic_excess_speed(self) -> float:
        return excess_speed(self.eccentricity, self.mu, self.semimajor_axis)

    @classmethod
    def from_orbital_elements(cls, e: float, rp: float, nu: float, mu: float,
                              raan: float = 0.0, i: float = 0.0, omega: float = 0.0):
        """
        Alternative constructor for the Orbit class. Takes the orbital elements and the
        gravitational parameter as arguments.

        Parameters
        ----------
        e : float
            Eccentricity
        rp : float
            Radius of periapsis (km)
        nu : float
            The true anomaly of the satellite (rad)
        mu : float
            The gravitational parameter of the central body (km^3/s^2)
        raan : float
            The right ascension of the ascending node (rad), default = 0.0
        i : float
            The orbital inclination (rad), default = 0.0
        omega : float
            The argument of periapsis (rad), default = 0.0

        Returns
        -------
        Orbit
            A new instance of Orbit
        """
        asymp_anomaly = asymptote_anomaly(e)
        nu_check = abs(nu)
        if np.isclose(nu_check, asymp_anomaly) or nu_check > asymp_anomaly:
            raise ValueError("State isn't defined at infinity")

        r, v = state_pf_from_e_rp(e, rp, nu, mu)

        eci_trans = perifocal_to_eci_trans_mat(raan, i, omega)
        r, v = np.matmul(eci_trans, r), np.matmul(eci_trans, v)

        return cls(r, v, mu)


