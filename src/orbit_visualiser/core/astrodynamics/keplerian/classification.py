from orbit_visualiser.core.astrodynamics.types import OrbitType, OrbitMotion
from math import pi

ECC_TOL = 1e-8

def orbit_type(e: float) -> OrbitType:
    """
    Returns the orbit type based on the eccentricity.

    Parameters
    ----------
    e : float
        Eccentricity

    Returns
    -------
    OrbitType
        The OrbitType enum
    """
    if e <= ECC_TOL:
        return OrbitType.CIRCULAR

    elif abs(e - 1.0) <= ECC_TOL :
        return OrbitType.PARABOLIC

    elif e < 1.0:
        return OrbitType.ELLIPTICAL

    return OrbitType.HYPERBOLIC

def orbit_motion_type(i: float) -> OrbitMotion:
    """
    Returns the orbit motion type based on the inclination.

    Parameters
    ----------
    i : float
        Inclination (rad)

    Returns
    -------
    OrbitMotion
        The OrbitMotion enum
    """
    if pi/2 < i <= pi:
        return OrbitMotion.RETROGRADE

    return OrbitMotion.PROGRADE