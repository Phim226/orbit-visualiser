import numpy as np
from orbit_visualiser.core.astrodynamics.types import OrbitType

def specific_ang_momentum(mu: float, p: float) -> float:
    """
    Calculates the specific angular momentum of a satellite using the gravitational parameter and
    the semi-parameter.

    Parameters
    ----------
    mu : float
        Gravitational parameter (km^3/s^2)
    p : float
        Semi-parameter (km)

    Returns
    -------
    float
        The specific angular momentum (km^2/s)
    """
    return np.sqrt(mu*p)

def specific_orbital_energy(orbit_type: OrbitType, mu: float, a: float) -> float:
    """
    Calculates the specific orbital energy using the gravitational parameter and semi-major axis.

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
        The mechanical energy per unit mass (km^2/s^2)
    """
    if orbit_type is OrbitType.PARABOLIC:
        return 0.0

    return -mu/(2*a)

def characteristic_energy(orbit_type: OrbitType, mu: float, a: float) -> float:
    """
    Calculates the characteristic energy from the semi-major axis and gravitational parameter.

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
        The characteristic energy (km^2/s^2)
    """
    if orbit_type is OrbitType.PARABOLIC:
        return 0.0

    c3 = mu/a
    if orbit_type in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
        c3 *= -1

    return c3

def excess_velocity(orbit_type: OrbitType, mu: float, a: float) -> float:
    """
    Calculates the hyperbolic excess velocity (the velocity magnitude at infinity) for open orbits
    from the semimajor axis and gravitational parameter.

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
        The hyperbolic excess velocity (km/s)
    """
    if orbit_type in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
        return np.nan

    return np.sqrt(mu/a)