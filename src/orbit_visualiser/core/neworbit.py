from dataclasses import dataclass
from typing import Sequence
import numpy as np
from numpy.typing import NDArray
from orbit_visualiser.core.astrodynamics.keplerian.state import state_pf_from_e_rp
from orbit_visualiser.core.astrodynamics.keplerian.elements import (eccentricity_from_state, semi_parameter_from_momentum,
                                                                    periapsis, semimajor_axis, semiminor_axis,
                                                                    apoapsis, asymptote_anomaly, turning_angle,
                                                                    aiming_radius, orbital_period, mean_motion,
                                                                    true_anomaly_from_state)
from orbit_visualiser.core.astrodynamics.keplerian.dynamics import (specific_orbital_energy, characteristic_energy,
                                                                    excess_velocity, specific_ang_momentum_from_state)
from orbit_visualiser.core.astrodynamics.keplerian.classification import orbit_type

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
    pos : Sequence | NDArray[np.float64]
        The perifocal position vector of the satellite (km)
    vel : Sequence | NDArray[np.float64]
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
            pos: Sequence | NDArray[np.float64],
            vel: Sequence | NDArray[np.float64],
            mu: float
    ):
        self._pos = pos
        self._vel = vel
        self._mu = mu

        self.compute_orbital_elements(pos, vel, mu)

    @property
    def position(self) -> Sequence | NDArray[np.float64]:
        return self._pos

    @property
    def velocity(self) -> Sequence | NDArray[np.float64]:
        return self._vel

    @property
    def gravitational_parameter(self) -> float:
        return self._mu


    @classmethod
    def from_orbital_elements(cls, e: float, rp: float, mu: float, nu: float):
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
        r, v = state_pf_from_e_rp(e, rp, mu, nu)
        return cls(r, v, mu)

    def compute_orbital_elements(
            self,
            pos: Sequence | NDArray[np.float64],
            vel: Sequence | NDArray[np.float64],
            mu: float
    ) -> None:
        """
        Computes the orbital elements, including geometric and energy invariants.

        Parameters
        ----------
        pos : Sequence | NDArray[np.float64]
            The perifocal position vector of the satellite (km)
        vel : Sequence | NDArray[np.float64]
            The perifocal velocity vector of the satellite (km/s)
        mu : float
            The gravitational parameter of the central body (km^3/s^2)
        """
        self.eccentricity = eccentricity_from_state(pos, vel, mu)
        self.true_anomaly = true_anomaly_from_state(pos)
        self.orbit_type = orbit_type(self.eccentricity)
        self.semi_parameter = semi_parameter_from_momentum(specific_ang_momentum_from_state(pos, vel), mu)
        self.radius_of_periapsis = periapsis(self.semi_parameter, self.eccentricity)
        self.semimajor_axis = semimajor_axis(self.orbit_type, self.eccentricity, self.radius_of_periapsis)
        self.semiminor_axis = semiminor_axis(self.orbit_type, self.eccentricity, self.semimajor_axis)
        self.radius_of_apoapsis = apoapsis(self.orbit_type, self.eccentricity, self.semimajor_axis)
        self.asymptote_anomaly = asymptote_anomaly(self.orbit_type, self.eccentricity)
        self.turning_angle = turning_angle(self.orbit_type, self.eccentricity)
        self.aiming_radius = aiming_radius(self.orbit_type, self.semiminor_axis)
        self.orbital_period = orbital_period(self.orbit_type, mu, self.semimajor_axis)
        self.mean_motion = mean_motion(self.orbit_type, self.orbital_period, mu, self.semi_parameter, self.semimajor_axis)
        self.specific_energy = specific_orbital_energy(self.orbit_type, mu, self.semimajor_axis)
        self.characteristic_energy = characteristic_energy(self.orbit_type, mu, self.semimajor_axis)
        self.hyperbolic_excess_velocity = excess_velocity(self.orbit_type, mu, self.semimajor_axis)

