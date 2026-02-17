from math import pi
import numpy as np
from numpy.typing import NDArray
from typing import Callable, Literal
from orbit_visualiser.core.astrodynamics.types import OrbitType
from orbit_visualiser.core.astrodynamics.keplerian.elements import semi_parameter_from_eccentricity
from orbit_visualiser.core.astrodynamics.keplerian.dynamics import specific_ang_momentum
from orbit_visualiser.core.astrodynamics.keplerian.classification import orbit_type

def state_pf_from_e_rp(
        e: float,
        rp: float,
        mu: float,
        nu: float,
        state: Literal["pos", "vel", "both"] = "both"
    ) -> list[NDArray[np.float64]]:
    """
    Takes the eccentricity, radius of periapsis, gravitational parameter and true anomaly, and returns
    the perifocal position and velocity vectors of a satellite at the given true anomaly on an
    analytical Keplerian orbit defined by the eccentricity and radius of periapsis.

    Parameters
    ----------
    e : float
        Eccentricity
    rp : float
        Radius of periapsis (km)
    mu : float
        Gravitational parameter (km^3/s^2)
    nu : float
        True anomaly (rads)
    state : Literal["pos", "vel", "both"]
        String literal indicating the state vector(s) being returned, default = both
    Returns
    -------
    list[NDArray[np.float64]]
        A list containing the numpy arrays of the perifocal position (km) and perifocal velocity (km/s)
    """
    p = semi_parameter_from_eccentricity(e, rp)
    h = specific_ang_momentum(mu, p)

    r = perifocal_position(e, p, nu)
    v = perifocal_velocity(e, mu, h, nu)

    if state == "pos":
        return [r]
    elif state == "vel":
        return [v]
    elif state == "both":
        return [r, v]

def perifocal_position(e: float, p: float, nu: float) -> NDArray[np.float64]:
    """
    Perifocal orbit equation.

    Parameters
    ----------
    e : float
        Eccentricity
    p : float
        Semi-parameter (km)
    nu : float
        True anomaly (rads)

    Returns
    -------
    NDArray[np.float64]
        The numpy array of the perifocal orbital position
    """
    return p*(1/(1 + e*np.cos(nu)))*np.array([np.cos(nu), np.sin(nu)])

def perifocal_velocity(e: float, mu: float, h: float, nu: float) -> NDArray[np.float64]:
    """
    The perifocal velocity equation.

    Parameters
    ----------
    e : float
        Eccentricity
    mu : float
        Gravitational parameter (km^3/s^2)
    h : float
        Specific angular momentum (km^2/s)
    nu : float
        True anomaly (rads)

    Returns
    -------
    NDArray[np.float64]
        The numpy array of the perifocal velocity
    """
    return (mu/h)*np.array([-np.sin(nu), e + np.cos(nu)])


def radial_azimuthal_velocity(
        nu: float,
        mu: float,
        h: float,
        e: float,
        asymp_anomaly: float
) -> NDArray[np.float64]:
    """
    Calculates the radial-azimuthal velocity vector using the gravitational parameter, angular
    momentum, eccentricity and true anomaly.

    Parameters
    ----------
    nu : float
        True anomaly (rads)
    mu : float
        Gravitational parameter (km^3/s^2)
    h : float
        Specific angular momentum (km^2/s)
    e : float
        Eccentricity
    asymp_anomaly : float
        The asymptote of the free anomaly (rads)

    Returns
    -------
    NDArray[np.float64]
        The radial-azimuthal velocity vector [v_radial, v_azimuthal] (km/s)
    """
    v_rad = 0.0 if orbit_type(e) is OrbitType.CIRCULAR else (mu/h)*e*np.sin(nu)
    v_azim = (0.0 if np.isclose(abs(nu), asymp_anomaly)
              else (mu/h)*(1 + e*np.cos(nu)))

    return np.array([v_rad, v_azim])

def speed(velocity: NDArray[np.float64]) -> float:
    """
    Velocity magnitude helper function.

    Parameters
    ----------
    velocity : NDArray[np.float64]
        Velocity vector (km/s)

    Returns
    -------
    float
        The magnitude of the velocity (km/s)
    """
    return np.linalg.norm(velocity)

def radius_from_state(r: NDArray[np.float64]) -> float:
    """
    Radial magnitude helper function.

    Parameters
    ----------
    r : NDArray[np.float64]
        Position vector (km)

    Returns
    -------
    float
        Magnitude of the radius vector (km)
    """
    return np.linalg.norm(r)

def radius_from_orbit_eq(nu: float, asymp_anomaly: float, p: float, e: float) -> float:
    """
    Calculates the magnitude of the radius using the orbit equation and the semi-parameter.

    Parameters
    ----------
    nu : float
        True anomaly (rads)
    asymp_anomaly : float
        The asymptote of the free anomaly (rads)
    p : float
        Semi-parameter (km)
    e : float
        Eccentricity

    Returns
    -------
    float
        The length of the radius vector (km)
    """
    if np.isclose(abs(nu), asymp_anomaly):
        return np.nan

    return p/(1 + e*np.cos(nu))

def escape_velocity(nu: float, asymp_anomaly: float, mu: float, r: float) -> float:
    """
    Calculates escape velocity at the given orbital radius and gravitational parameter.

    Parameters
    ----------
    nu : float
        The true anomaly (rads)
    asymp_anomaly : float
        The true anomaly of the asymptot (rads)
    mu : float
        Gravitational parameter (km^3/s^2)
    r : float
        Orbital radius (km)

    Returns
    -------
    float
        The escape velocity (km/s)
    """
    if np.isclose(abs(nu), asymp_anomaly):
        return 0.0

    return np.sqrt(2*mu/r)

def flight_angle(nu: float, asymp_anomaly: float, e: float) -> float:
    """
    Calculates the flight angle of a satellite at the given true anomaly and eccentricity.

    Parameters
    ----------
    nu : float
        True anomaly (rads)
    asymp_anomaly : float
        The true anomaly of the asymptot (rads)
    e : float
        Eccentricity

    Returns
    -------
    float
        The satellite flight angle (rads)
    """
    if np.isclose(abs(nu), asymp_anomaly):
        return np.sign(nu)*pi/2

    return np.arctan2(e*np.sin(nu), 1 + e*np.cos(nu))

def time_since_periapsis(m_anomaly: float, period: float, p: float, h: float, e: float) -> float:
    """
    Calculates the time since periapsis, using different formulae based on orbit type.

    Parameters
    ----------
    m_anomaly : float
        The mean anomaly (rads)
    period : float
        The orbital period (s)
    p : float
        Semi-parameter (km)
    h : float
        Specific angular momentum (km^2/s)
    e : float
        Eccentricity

    Returns
    -------
    float
        Time since periapsis (s)
    """
    type = orbit_type(e)
    if type in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
        return period*m_anomaly/(2*pi)

    elif type is OrbitType.PARABOLIC:
        return (p**2/h)*m_anomaly

    elif type is OrbitType.HYPERBOLIC:
        return (p**2/h)*(e**2 - 1)**(-1.5)*m_anomaly