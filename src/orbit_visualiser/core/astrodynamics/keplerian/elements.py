from math import pi
import numpy as np
from numpy.typing import NDArray
from orbit_visualiser.core.astrodynamics.types import OrbitType
from orbit_visualiser.core.astrodynamics.keplerian.classification import orbit_type

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
        The eccentricity vector
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
    r_norm = np.linalg.norm(r)
    v_norm = np.linalg.norm(v)
    return np.sqrt((v_norm**2*r_norm - mu)**2 + (2*mu - v_norm**2*r_norm)*np.dot(r, v)**2/r_norm)/mu

def true_anomaly(r: NDArray[np.float64], e: NDArray[np.float64], v_r: float) -> float:
    """
    Calculates the true anomaly from the current position vector, the eccentricity vector and the
    radial speed.

    Parameters
    ----------
    r : NDArray[np.float64]
        Position vector of the satellite (km)
    e : NDArray[np.float64]
        The eccentricity vector
    v_r : float
        The radial speed (km/s)
    Returns
    -------
    float
        True anomaly (rads)
    """
    per_vect = e
    per_vect_norm = np.linalg.norm(e)

    # per_vect stands for periapsis vector. In the case where e = 0 then we set the periapsis vector
    # to be the x axis, so that the argument of periapsis is 0. (Will break once inclination is considered properly)
    if np.isclose(per_vect_norm, 0):
        per_vect = np.array([1.0, 0.0, 0.0])
        per_vect_norm = np.linalg.norm(per_vect)

    true_anomaly = np.arccos(np.dot(per_vect, r)/(per_vect_norm*np.linalg.norm(r)))

    if (v_r < 0 and not np.isclose(v_r, 0)) or (np.isclose(np.linalg.norm(e), 0) and r[1] < 0):
        true_anomaly = 2*pi - true_anomaly

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

def semimajor_axis(e: float, rp: float) -> float:
    """
    Calculates the semimajor axis using the eccentricity and radius of periapsis.

    Parameters
    ----------
    e : float
        Eccentricity
    rp : float
        Radius of periapsis (km)

    Returns
    -------
    float
        Semi-major axis (km)
    """
    if orbit_type(e) is OrbitType.PARABOLIC:
        return np.nan

    return rp/(1 - e)

def semiminor_axis(e: float, a: float) -> float:
    """
    Calculates the semi-minor axis using the eccentricity and semi-major axis.

    Parameters
    ----------
    e : float
        Eccentricity
    a : float
        Semi-major axis (km)

    Returns
    -------
    float
        Semi-minor axis (km)
    """
    type = orbit_type(e)
    if type in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
        return a*np.sqrt(1 - e**2)

    elif type is OrbitType.HYPERBOLIC:
        return a*np.sqrt(e**2 - 1)

    return np.nan

def radius_of_periapsis(p: float, e: float) -> float:
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

def periapsis(p: float, e: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    Calculate the periapsis vector from the eccentricity vector and the semi-parameter.

    Parameters
    ----------
    p : float
        The semi-parameter (km)
    e : NDArray[np.float64]
        The eccentricity vector

    Returns
    -------
    NDArray[np.float64]
        Periapsis vector (km)
    """
    e_norm = np.linalg.norm(e)
    return p/((e_norm(1 + e_norm)))*e

def radius_of_apoapsis(e: float, a: float) -> float:
    """
    Calculates the radius of apoapsis using the eccentricity and semi-major axis.

    Parameters
    ----------
    e : float
        Eccentricity
    a : float
        Semi-major axis (km)

    Returns
    -------
    float
        Radius of apoapsis (km)
    """
    if orbit_type(e) is OrbitType.PARABOLIC:
        return np.nan

    return a*(1 + e)

def apoapsis(e: NDArray[np.float64], a: float) -> NDArray[np.float64]:
    """
    Calculate the apoapsis vector from the eccentricity vector and the semi-major axis.

    Parameters
    ----------
    e : NDArray[np.float64]
        Eccentricity vector
    a : float
        Semi-major axis (km)

    Returns
    -------
    NDArray[np.float64]
        Apoapsis vector (km)
    """

    e_norm = np.linalg.norm(e)
    if orbit_type(e_norm) is OrbitType.PARABOLIC:
        return np.empty(e.shape)

    return -(a*(1 + e_norm)/e_norm)*e

def asymptote_anomaly(e: float) -> float:
    """
    Calculate the asymptote of the true anomaly for open orbits using the eccentricity.

    Parameters
    ----------
    e : float
        Eccentricity

    Returns
    -------
    float
        The true anomaly of the asymptote (rads)
    """
    if orbit_type(e) in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
        return np.nan

    e = 1 if np.isclose(e, 1) else e

    return np.arccos(-1/e)

def turning_angle(e: float) -> float:
    """
    Calculate the turning angle for open orbits using the eccentricity

    Parameters
    ----------
    e : float
        Eccentricity

    Returns
    -------
    float
        The turning angle (rads)
    """
    if orbit_type(e) in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
        return np.nan

    e = 1 if np.isclose(e, 1) else e

    return 2*np.arcsin(1/e)

def aiming_radius(b: float) -> float:
    """
    Calculate the aiming radius for hyperbolic orbits using the semi-minor axis

    Parameters
    ----------
    b : float
        Semi-minor axis (km)

    Returns
    -------
    float
        Aiming radius (km)
    """
    if b > 0 or np.isclose(b, np.nan, equal_nan = True):
        return np.nan

    return -b

def orbital_period(e: float, mu: float, a: float) -> float:
    """
    Calculates the orbital period of closed orbits.

    Parameters
    ----------
    e : float
        The eccentricity
    mu : float
        Gravitational parameter (km^3/s^2)
    a : float
        Semi-major axis (km)

    Returns
    -------
    float
        The orbital period (s)
    """
    if orbit_type(e) not in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
        return np.nan

    return (2*pi/np.sqrt(mu))*np.sqrt(a)**3

def mean_motion(e: float, mu: float, p: float, a: float) -> float:
    """
    Calculates the mean motion, using different formulae based on orbit type.

    Parameters
    ----------
    e : float
        The eccentricity
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
    if orbit_type(e) is OrbitType.PARABOLIC:
        return 2*np.sqrt(mu/(p**3))

    return np.sqrt(mu/abs(a**3))

def inclination(h: NDArray[np.float64]) -> float:
    """
    Calculates the orbital inclination from the specific angular momentum vector.

    Parameters
    ----------
    h : NDArray[np.float64]
        Specific angular momentum (km^2/s)

    Returns
    -------
    float
        The orbital inclination (rad)
    """
    return np.arccos(h[2]/np.linalg.norm(h))

def node_line(h: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    Calculates the vector defining the node line from the specific angular momentum.

    Parameters
    ----------
    h : NDArray[np.float64]
        Specific angular momentum (km^2/s)

    Returns
    -------
    NDArray[np.float64]
        The node line vector
    """
    node_line = np.cross([0, 0, 1], h)
    if np.allclose(node_line, np.zeros(node_line.shape)):
        node_line = np.array([1.0, 0.0, 0.0])

    return node_line

def right_ascen_of_ascending_node(node_line: NDArray[np.float64]) -> float:
    """
    Calculates the right ascension of the ascending node from the node line vector.
    Returns nan if node_line is 0.

    Parameters
    ----------
    node_line : NDArray[np.float64]
        Node line vector

    Returns
    -------
    float
        Right ascension of the ascending node (rad)
    """
    raan = np.arccos(node_line[0]/np.linalg.norm(node_line))

    if node_line[1] < 0:
        raan = 2*pi - raan

    return raan

def argument_of_periapsis(node_line: NDArray[np.float64], e: NDArray[np.float64]) -> float:
    """
    Calculate the argument of periapsis from the node line vector and eccentricity vector.
    Returns nan if node_line is 0.

    Parameters
    ----------
    node_line : NDArray[np.float64]
        The node line vector
    e : NDArray[np.float64]
        The eccentricity vector

    Returns
    -------
    float
        The argument of periapsis (rad)
    """
    if np.allclose(e, np.zeros(e.shape)):
        e = np.array([1.0, 0.0, 0.0])

    arg_periapsis = np.arccos(np.dot(node_line, e)/(np.linalg.norm(node_line)*np.linalg.norm(e)))

    if e[2] < 0:
        arg_periapsis = 2*pi - arg_periapsis

    return arg_periapsis