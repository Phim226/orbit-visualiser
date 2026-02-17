from orbit_visualiser.core.astrodynamics.types import OrbitType

ECC_TOL = 1e-8

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
    if e <= ECC_TOL:
        return OrbitType.CIRCULAR

    elif abs(e - 1.0) <= ECC_TOL :
        return OrbitType.PARABOLIC

    elif e < 1.0:
        return OrbitType.ELLIPTICAL

    return OrbitType.HYPERBOLIC