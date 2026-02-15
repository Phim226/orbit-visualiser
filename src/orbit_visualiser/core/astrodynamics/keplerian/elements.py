from math import pi
import numpy as np
from numpy.typing import NDArray
from orbit_visualiser.core.astrodynamics.types import OrbitType

def eccentricity_vector_from_state(r: NDArray[np.float64], v: NDArray[np.float64], mu: float) -> NDArray[np.float64]:
    """
    Calculates the eccentricity vector from the position and velocity vectors of a satellite and
    the gravitational parameter.

    Parameters
    ----------
    r : NDArray[np.float64]
        Position vector of the satellite (km)
    v : NDArray[np.float64]
        Velocity vector of the satellite (km/s)
    mu : float
        Gravitational parameter (km^3/s^2)

    Returns
    -------
    NDArray[np.float64]
        _description_
    """
    return (1/mu)*((np.linalg.norm(v)**2 - mu/np.linalg.norm(r))*r - np.dot(r, v)*v)

def eccentricity_from_state(r: NDArray[np.float64], v: NDArray[np.float64], mu: float) -> float:
    """
    Calculates the eccentricity from the position and velocity vectors of a satellite and
    the gravitational parameter.

    Parameters
    ----------
    r : NDArray[np.float64]
        Position vector of the satellite (km)
    v : NDArray[np.float64]
        Velocity vector of the satellite (km/s)
    mu : float
        Gravitational parameter (km^3/s^2)

    Returns
    -------
    float
        Eccentricity
    """
    return np.linalg.norm(eccentricity_vector_from_state(r, v, mu))

def true_anomaly_from_state(r: NDArray[np.float64]) -> float:
    """
    Calculates the true anomaly from the current position vector.

    Parameters
    ----------
    r : NDArray[np.float64]
        Position vector of the satellite (km)

    Returns
    -------
    float
        True anomaly (rads)
    """
    true_anomaly = np.atan2(r[1], r[0])
    if true_anomaly < 0:
        return true_anomaly + 2*pi

    return true_anomaly

def semi_parameter_from_momentum(h: float, mu: float) -> float:
    """
    Calculates the semi-parameter using the gravitational parameter and specific angular momentum.

    Parameters
    ----------
    h : float
        Specific angular momentum (km^2/s)
    mu : float
        Gravitational parameter (km^3/s^2)

    Returns
    -------
    float
        Semi-parameter (km)
    """
    return h**2/mu

def semi_parameter_from_eccentricity(e: float, rp: float) -> float:
    """
    Calculates the semi-parameter using the eccentricity and radius of periapsis.

    Parameters
    ----------
    e : float
        Eccentricity
    rp : float
        Radius of periapsis (km)

    Returns
    -------
    float
        Semi-parameter (km)
    """
    return rp*(1 + e)

def semimajor_axis(orbit_type: OrbitType, e: float, rp: float) -> float:
    """
    Calculates the semimajor axis using the eccentricity and radius of periapsis.

    Parameters
    ----------
    orbit_type : OrbitType
        The orbit type enum
    e : float
        Eccentricity
    rp : float
        Radius of periapsis (km)

    Returns
    -------
    float
        Semi-major axis (km)
    """
    if orbit_type is OrbitType.PARABOLIC:
        return np.inf

    return rp/(1 - e)

def semiminor_axis(orbit_type: OrbitType, e: float, a: float) -> float:
    """
    Calculates the semi-minor axis using the eccentricity and semi-major axis.

    Parameters
    ----------
    orbit_type : OrbitType
        The orbit type enum
    e : float
        Eccentricity
    a : float
        Semi-major axis (km)

    Returns
    -------
    float
        Semi-minor axis (km)
    """
    if orbit_type in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
        return a*np.sqrt(1 - e**2)

    elif orbit_type is OrbitType.HYPERBOLIC:
        return a*np.sqrt(e**2 - 1)

    return np.inf

def periapsis(p: float, e: float) -> float:
    """
    Calculates the radius of periapsis from the semi-parameter and eccentricity.

    Parameters
    ----------
    p : float
        Semi-parameter (km)
    e : float
        Eccentricity

    Returns
    -------
    float
        Radius of periapsis (km)
    """
    return p/(1 + e)

def apoapsis(orbit_type: OrbitType, e: float, a: float) -> float:
    """
    Calculates the radius of apoapsis using the eccentricity and semi-major axis.

    Parameters
    ----------
    orbit_type : OrbitType
        The orbit type enum
    e : float
        Eccentricity
    a : float
        Semi-major axis (km)

    Returns
    -------
    float
        Radius of apoapsis (km)
    """
    if orbit_type is OrbitType.PARABOLIC:
        return np.inf

    return a*(1 + e)

def asymptote_anomaly(orbit_type: OrbitType, e: float) -> float:
    """
    Calculate the asymptote of the true anomaly for open orbits using the eccentricity.

    Parameters
    ----------
    orbit_type : OrbitType
        The orbit type enum
    e : float
        Eccentricity

    Returns
    -------
    float
        The true anomaly of the asymptote (rads)
    """
    if orbit_type in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
        return np.nan

    return np.arccos(-1/e)

def turning_angle(orbit_type: OrbitType, e: float) -> float:
    """
    Calculate the turning angle for open orbits using the eccentricity

    Parameters
    ----------
    orbit_type : OrbitType
        The orbit type enum
    e : float
        Eccentricity

    Returns
    -------
    float
        The turning angle (rads)
    """
    if orbit_type in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
        return np.nan

    return 2*np.arcsin(1/e)

def aiming_radius(orbit_type: OrbitType, b: float) -> float:
    """
    Calculate the aiming radius for hyperbolic orbits using the semi-minor axis

    Parameters
    ----------
    orbit_type : OrbitType
        The orbit type enum
    b : float
        Semi-minor axis (km)

    Returns
    -------
    float
        Aiming radius (km)
    """
    if orbit_type is not OrbitType.HYPERBOLIC:
        return np.nan

    return -b

def orbital_period(orbit_type: OrbitType, mu: float, a: float) -> float:
    """
    Calculates the orbital period of closed orbits.

    Parameters
    ----------
    orbit_type : OrbitType
        The orbit type enum
    mu : float
        Gravitational parameter (km^3/s^2)
    a : float
        Semi-major axis (km)

    Returns
    -------
    float
        The orbital period (s)
    """
    if orbit_type not in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
        return np.nan

    return (2*pi/np.sqrt(mu))*np.sqrt(a)**3

def mean_motion(orbit_type: OrbitType, period: float, mu: float, p: float, a: float) -> float:
    """
    Calculates the mean motion, using different formulae based on orbit type.

    Parameters
    ----------
    orbit_type : OrbitType
        The orbit type enum
    period : float
        The orbital period (s)
    mu : float
        Gravitational parameter (km^3/s^2)
    p : float
        Semi-parameter (km)
    a : float
        Semi-major axis (km)

    Returns
    -------
    float
        The mean motion (rads/s)
    """
    if orbit_type in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL) :
        return 2*pi/period

    elif orbit_type is OrbitType.PARABOLIC:
        return 2*np.sqrt(mu/p)

    elif orbit_type is OrbitType.HYPERBOLIC:
        return 2*np.sqrt(mu/abs(a**3))