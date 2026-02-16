import numpy as np
from numpy.typing import NDArray
from typing import Callable, Sequence
from math import pi
from orbit_visualiser.core.astrodynamics.keplerian.state import (perifocal_position_eq, perifocal_velocity_eq,
                                                                 speed, radial_azimuthal_velocity, escape_velocity,
                                                                 radius_from_state, flight_angle, time_since_periapsis)
from orbit_visualiser.core.astrodynamics.keplerian.anomlies import mean_anomaly, eccentric_anomaly
from orbit_visualiser.core.astrodynamics.keplerian.dynamics import specific_ang_momentum_from_state
from orbit_visualiser.core.astrodynamics.keplerian.elements import true_anomaly_from_state
from orbit_visualiser.core.orbit import Orbit, CentralBody
from orbit_visualiser.core.neworbit import NewOrbit

class NewSatellite():
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
    def orbit(self) -> NewOrbit:
        return NewOrbit(self._pos, self._vel, self._central_body.mu)

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
            orbit.orbit_type,
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
            orbit.orbit_type,
            orbit.eccentricity,
            self.true_anomaly
        )

    @property
    def mean_anomaly(self) -> float:
        orbit = self.orbit
        return mean_anomaly(
            orbit.orbit_type,
            orbit.eccentricity,
            self.true_anomaly,
            self.eccentric_anomaly
        )

    @property
    def time_since_periapsis(self) -> float:
        orbit = self.orbit
        return time_since_periapsis(
            orbit.orbit_type,
            self.mean_anomaly,
            orbit.orbital_period,
            orbit.semi_parameter,
            self.specific_angular_momentum,
            orbit.eccentricity
        )

# TODO: Split formulae from Satellite class.
# TODO: Update satellite properties when true anomaly is set rather than calling it explicitly outside the class
class Satellite():


    def __init__(self, orbit: Orbit, central_body: CentralBody, nu: float = 0.00):
        self._orbit = orbit
        self._central_body = central_body
        self._nu = nu # true anomaly in rads

        self.update_satellite_properties()

    @property
    def nu(self) -> float:
        return self._nu

    @nu.setter
    def nu(self, nu: float) -> None:
        self._nu = nu

    @property
    def pos_pf_eq(self) -> Callable[[float], NDArray[np.float64]]:
        return self._pos_pf_eq

    @property
    def vel_pf_eq(self) -> Callable[[float], NDArray[np.float64]]:
        return self._vel_pf_eq

    @property
    def pos_pf(self) -> NDArray:
        return self._pos_pf

    @property
    def vel_pf(self) -> NDArray:
        return self._vel_pf

    @property
    def h(self) -> float:
        return self._h

    @property
    def v_azim(self) -> float:
        return self._v_azim

    @property
    def v_radial(self) -> float:
        return self._v_radial

    @property
    def v(self) -> float:
        return self._v

    @property
    def v_esc(self) -> float:
        return self._v_esc

    @property
    def v_inf(self) -> float:
        return self._v_inf

    @property
    def r(self) -> float:
        return self._r

    @property
    def eps(self) -> float:
        return self._eps

    @property
    def gam(self) -> float:
        return self._gam

    @property
    def c3(self) -> float:
        return self._c3

    @property
    def period(self) -> float:
        return self._period

    @property
    def n(self) -> float:
        return self._n

    @property
    def e_anomaly(self) -> float:
        return self._e_anomaly

    @property
    def m_anomaly(self) -> float:
        return self._m_anomaly

    @property
    def t_p(self) -> float:
        return self._t_p

    def update_satellite_properties(self) -> None:
        nu, mu, e, p = self._nu, self._central_body.mu, self._orbit.e, self._orbit.p
        h, t_asymp, a = self._specific_ang_momentum(mu, p), self._orbit.t_asymp, self._orbit.a
        self._h = h

        # Perifocal equations
        self._pos_pf_eq = perifocal_position_eq(e, p)
        self._vel_pf_eq = perifocal_velocity_eq(e, mu, h)

        # Helper quantities
        den = 1 + e*np.cos(nu) # Denominator of the orbit equation
        mu_over_h = mu/h

        # Geometry
        self._pos_pf = self._pf_position(nu, t_asymp)
        self._r = self._radius(nu, t_asymp, p, den)

        # Kinematics
        vel_pf = self._pf_velocity(nu)
        self._vel_pf = vel_pf
        self._v_azim = self._azimuthal_velocity(nu, t_asymp, mu_over_h, den)
        self._v_radial = self._radial_velocity(mu_over_h, e, nu)
        self._v = self._velocity(vel_pf[0], vel_pf[1])
        self._v_esc = self._escape_velocity(nu, mu, self._r, t_asymp)
        self._v_inf = self._excess_velocity(mu, abs(a))
        self._gam = self._flight_angle(e, nu, t_asymp, den)

        # Energy
        self._eps = self._specific_energy(mu, a)
        self._c3 = self._characteristic_energy(mu, abs(a))

        # Time
        period = self._orbital_period(mu, a)
        self._period = period
        self._n = self._mean_motion(period, mu, p, a)
        e_anomaly = self._eccentric_anomaly(e, nu, t_asymp)
        self._e_anomaly = e_anomaly
        m_anomaly = self._mean_anomaly(e, nu, e_anomaly, t_asymp)
        self._m_anomaly = m_anomaly
        self._t_p = self._time_since_periapsis(m_anomaly, period, p, h, e)

    @staticmethod
    def _specific_ang_momentum(mu: float, p: float) -> float:
        return np.sqrt(mu*p)

    def _pf_position(self, nu: float, t_asymp: float) -> NDArray[np.float64]:
        pos_eq = self._pos_pf_eq
        return pos_eq(nu, t_asymp)

    def _pf_velocity(self, nu: float) -> NDArray[np.float64]:
        return self._vel_pf_eq(nu)

    def _azimuthal_velocity(self, nu: float, t_asymp: float, mu_over_h: float, den: float) -> float:
        if np.isclose(abs(nu), t_asymp, atol = 0.0001, rtol = 0):
            return 0.0

        return mu_over_h*den

    def _radial_velocity(self, mu_over_h: float, e: float, nu: float) -> float:
        if self._orbit.orbit_type == "circular":
            return 0.0

        return mu_over_h*e*np.sin(nu)

    def _velocity(self, v_x: float, v_y: float) -> float:
        return np.hypot(v_x, v_y)

    def _radius(self, nu: float, t_asymp: float, p: float, den: float) -> float:
        if np.isclose(abs(nu), t_asymp, atol = 0.0001, rtol = 0):
            return np.inf

        return p/den

    def _specific_energy(self, mu: float, a: float) -> float:
        if self._orbit.orbit_type == "parabolic":
            return 0.0

        return -mu/(2*a)

    def _escape_velocity(self, nu: float, mu: float, r: float, t_asymp: float) -> float:
        if np.isclose(abs(nu), t_asymp, atol = 0.0001, rtol = 0):
            return 0.0

        return np.sqrt(2*mu/r)

    def _excess_velocity(self, mu: float, a: float) -> float:
        if self._orbit.is_closed:
            return np.nan

        return np.sqrt(mu/a)

    def _flight_angle(self, e: float, nu: float, t_asymp: float, den: float) -> float:
        if np.isclose(abs(nu), t_asymp, atol = 0.0001, rtol = 0):
            return (nu/abs(nu))*pi/2

        return np.arctan2(e*np.sin(nu), den)

    def _characteristic_energy(self, mu: float, a: float) -> float:
        if self._orbit.orbit_type == "parabolic":
            return 0.0

        c3 = mu/a
        if self._orbit.is_closed:
            c3 *= -1

        return c3

    def _orbital_period(self, mu: float, a: float) -> float:
        if self._orbit.is_closed:
            return (2*pi/np.sqrt(mu))*np.sqrt(a)**3

        return np.nan

    def _mean_motion(self, period: float, mu: float, p: float, a: float) -> float:
        if self._orbit.is_closed:
            return 2*pi/period

        orbit_type = self._orbit.orbit_type
        if orbit_type == "parabolic":
            return 2*np.sqrt(mu/p)

        elif orbit_type == "hyperbolic":
            return 2*np.sqrt(mu/abs(a**3))

    def _eccentric_anomaly(self, e: float, nu: float, t_asymp: float):
        orbit_type = self._orbit.orbit_type
        if orbit_type == "circular":
            return nu

        elif orbit_type == "elliptical":
            e_anomaly = 2*np.arctan(np.sqrt((1 - e)/(1 + e))*np.tan(nu/2))
            if e_anomaly < 0:
                return e_anomaly + 2*pi

            return e_anomaly

        elif orbit_type == "parabolic":
            return np.nan

        elif orbit_type == "hyperbolic":
            if np.isclose(abs(nu), t_asymp, atol = 0.0001, rtol = 0):
                return (nu/abs(nu))*np.inf

            return 2*np.arctanh(np.sqrt((e - 1)/(e + 1))*np.tan(nu/2))


    def _mean_anomaly(self, e: float, nu: float, e_anomaly: float, t_asymp: float):
        if self._orbit.is_closed:
            return e_anomaly - e*np.sin(e_anomaly)

        orbit_type = self._orbit.orbit_type
        if np.isclose(abs(nu), t_asymp, atol = 0.0001, rtol = 0):
            return (nu/abs(nu))*np.inf

        elif orbit_type == "parabolic":
            return 0.5*np.tan(nu/2) + (1/6)*(np.tan(nu/2))**3

        elif orbit_type == "hyperbolic":
            return e*np.sinh(e_anomaly) - e_anomaly


    def _time_since_periapsis(self, m_anomaly: float, period: float, p: float, h: float, e: float) -> float:
        if np.isneginf(m_anomaly):
            return -np.inf

        elif np.isinf(m_anomaly):
            return np.inf

        if self._orbit.is_closed:
            return period*m_anomaly/(2*pi)

        orbit_type = self._orbit.orbit_type
        if orbit_type == "parabolic":
            return (p**2/h)*m_anomaly

        elif orbit_type == "hyperbolic":
            return (p**2/h)*(e**2 - 1)**(-1.5)*m_anomaly