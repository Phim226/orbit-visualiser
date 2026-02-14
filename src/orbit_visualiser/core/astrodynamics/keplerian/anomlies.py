import numpy as np
from math import pi
from orbit_visualiser.core.astrodynamics.types import OrbitType

def eccentric_anomaly(orbit_type: OrbitType, e: float, nu: float, asymp_anomaly: float) -> float:
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
    asymp_anomaly : float
        True anomaly of the asymptote (rads)

    Returns
    -------
    float
        The eccentric anomaly (rads)
    """
    if orbit_type is OrbitType.CIRCULAR:
        return nu

    elif orbit_type is OrbitType.ELLIPTICAL:
        e_anomaly = 2*np.arctan(np.sqrt((1 - e)/(1 + e))*np.tan(nu/2))
        if e_anomaly < 0:
            return e_anomaly + 2*pi

        return e_anomaly

    elif orbit_type is OrbitType.PARABOLIC:
        return np.nan

    elif orbit_type is OrbitType.HYPERBOLIC:
        if np.isclose(abs(nu), asymp_anomaly, atol = 0.0001, rtol = 0):
            return (nu/abs(nu))*np.inf

        return 2*np.arctanh(np.sqrt((e - 1)/(e + 1))*np.tan(nu/2))


def mean_anomaly(orbit_type: OrbitType, e: float, nu: float, e_anomaly: float, asymp_anomaly: float) -> float:
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
    asymp_anomaly : float
        True anomaly of the asymptote (rads)

    Returns
    -------
    float
        The mean anomaly (rads)
    """
    if np.isclose(abs(nu), asymp_anomaly, atol = 0.0001, rtol = 0):
        return (np.sign(nu))*np.inf

    if orbit_type in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
        return e_anomaly - e*np.sin(e_anomaly)

    elif orbit_type is OrbitType.PARABOLIC:
        return 0.5*np.tan(nu/2) + (1/6)*(np.tan(nu/2))**3

    elif orbit_type is OrbitType.HYPERBOLIC:
        return e*np.sinh(e_anomaly) - e_anomaly