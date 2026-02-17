import numpy as np
from numpy.typing import NDArray
from typing import Sequence
from orbit_visualiser.core.astrodynamics.keplerian.state import (speed, radial_azimuthal_velocity, escape_velocity,
                                                                 radius_from_state, flight_angle, time_since_periapsis)
from orbit_visualiser.core.astrodynamics.keplerian.anomlies import mean_anomaly, eccentric_anomaly
from orbit_visualiser.core.astrodynamics.keplerian.dynamics import specific_ang_momentum_from_state
from orbit_visualiser.core.astrodynamics.keplerian.elements import true_anomaly_from_state
from orbit_visualiser.core.orbit import CentralBody, Orbit

class Satellite():
    """
    Represents a satellite orbiting a central body. Tracks the dynamical state.

    Parameters
    ----------
    position : Sequence | NDArray[np.float64]
        The initial position vector of the satellite (km)
    velocity : float
        The initial velocity vector of the satellite (km/s)
    central_body: CentralBody
        The CentralBody object representing the body that the satellite is orbiting
    """

    def __init__(
            self,
            position: Sequence | NDArray[np.float64],
            velocity: Sequence | NDArray[np.float64],
            central_body: CentralBody
        ):
        if not np.all(np.isfinite(position)) or not np.all(np.isfinite(velocity)):
            raise ValueError("State vectors must only contain finite values.")

        self._pos = position
        self._vel = velocity
        self._central_body = central_body

    @property
    def position(self) -> Sequence | NDArray[np.float64]:
        """
        Perifocal position of the satellite.

        Returns
        -------
        Sequence | NDArray[np.float64]
            Perifocal position vector (km)
        """
        return self._pos

    @position.setter
    def position(self, new_pos: Sequence | NDArray[np.float64]) -> None:
        self._pos = new_pos

    @property
    def velocity(self) -> Sequence | NDArray[np.float64]:
        """
        Perifocal velocity of the satellite.

        Returns
        -------
        Sequence | NDArray[np.float64]
            Perifocal velocity vector (km/s)
        """
        return self._vel

    @velocity.setter
    def velocity(self, new_vel: Sequence | NDArray[np.float64]) -> None:
        self._vel = new_vel

    @property
    def central_body(self) -> CentralBody:
        return self._central_body

    @property
    def orbit(self) -> Orbit:
        return Orbit(self._pos, self._vel, self._central_body.mu)

    @property
    def radius(self) -> float:
        return radius_from_state(self._pos)

    @property
    def true_anomaly(self) -> float:
        return true_anomaly_from_state(self.position)

    @property
    def speed(self) -> float:
        return speed(self._vel)

    @property
    def radial_azimuthal_velocity(self) -> NDArray[np.float64]:
        orbit = self.orbit
        return radial_azimuthal_velocity(
            self.true_anomaly,
            self.central_body.mu,
            self.specific_angular_momentum,
            orbit.eccentricity,
            orbit.asymptote_anomaly
        )

    @property
    def flight_angle(self) -> float:
        orbit = self.orbit
        return flight_angle(self.true_anomaly, orbit.asymptote_anomaly, orbit.eccentricity)

    @property
    def escape_velocity(self) -> float:
        orbit = self.orbit
        return escape_velocity(
            self.true_anomaly,
            orbit.asymptote_anomaly,
            self.central_body.mu,
            self.radius
        )

    @property
    def specific_angular_momentum(self) -> float:
        return specific_ang_momentum_from_state(self._pos, self._vel)

    @property
    def eccentric_anomaly(self) -> float:
        orbit = self.orbit
        return eccentric_anomaly(
            orbit.eccentricity,
            self.true_anomaly
        )

    @property
    def mean_anomaly(self) -> float:
        orbit = self.orbit
        return mean_anomaly(
            orbit.eccentricity,
            self.true_anomaly,
            self.eccentric_anomaly
        )

    @property
    def time_since_periapsis(self) -> float:
        orbit = self.orbit
        return time_since_periapsis(
            self.mean_anomaly,
            orbit.orbital_period,
            orbit.semi_parameter,
            self.specific_angular_momentum,
            orbit.eccentricity
        )