from orbit_visualiser.core.astrodynamics.types import OrbitType

def orbit_type(e: float) -> OrbitType:
    """
    Returns the orbit type bases on the eccentricity.

    Parameters
    ----------
    e : float
        Eccentricity

    Returns
    -------
    OrbitType
        The OrbitType enum
    """
    if e == 0:
        return OrbitType.CIRCULAR

    elif 0 < e < 1:
        return OrbitType.ELLIPTICAL

    elif e == 1:
        return OrbitType.PARABOLIC

    elif e > 1:
        return OrbitType.HYPERBOLIC