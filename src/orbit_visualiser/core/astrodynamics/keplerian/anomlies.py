import numpy as np
from math import pi
from orbit_visualiser.core.astrodynamics.types import OrbitType

def eccentric_anomaly(orbit_type: OrbitType, e: float, nu: float) -> float:
    """
    Calculates the eccentric anomaly, using different formulae based on orbit type.

    Parameters
    ----------
    orbit_type : OrbitType
        The orbit type enum
    e : float
        Eccentricity
    nu : float
        True anomaly (rads)

    Returns
    -------
    float
        The eccentric anomaly (rads)
    """
    if orbit_type is OrbitType.CIRCULAR:
        return nu

    elif orbit_type is OrbitType.ELLIPTICAL:
        e_anomaly = 2*np.arctan(np.sqrt((1 - e)/(1 + e))*np.tan(nu/2))
        e_anomaly = np.where(e_anomaly < 0, e_anomaly + 2*pi, e_anomaly)

        return e_anomaly

    elif orbit_type is OrbitType.PARABOLIC:
        return np.nan

    elif orbit_type is OrbitType.HYPERBOLIC:
        return 2*np.arctanh(np.sqrt((e - 1)/(e + 1))*np.tan(nu/2))


def mean_anomaly(orbit_type: OrbitType, e: float, nu: float, e_anomaly: float) -> float:
    """
    Calculates the mean anomaly, using different formulae based on orbit type.

    Parameters
    ----------
    orbit_type : OrbitType
        The orbit type enum
    e : float
        Eccentricity
    nu : float
        True anomaly (rads)
    e_anomaly : float
        The eccentric anomaly (rads)


    Returns
    -------
    float
        The mean anomaly (rads)
    """
    if orbit_type in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
        return e_anomaly - e*np.sin(e_anomaly)

    elif orbit_type is OrbitType.PARABOLIC:
        return 0.5*np.tan(nu/2) + (1/6)*(np.tan(nu/2))**3

    elif orbit_type is OrbitType.HYPERBOLIC:
        return e*np.sinh(e_anomaly) - e_anomaly