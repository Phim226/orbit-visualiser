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

def orbit_closed(orbit_type: OrbitType) -> bool:
    """
    Returns a boolean for whether the orbit is closed. Circular and elliptical orbits are closed
    so return True, parabolic and hyperbolic are open so return False.

    Parameters
    ----------
    orbit_type : OrbitType
        The orbit type enum

    Returns
    -------
    bool
        Boolean value representing the closure of the orbit
    """
    closure_dict: dict[OrbitType, bool] = {
        OrbitType.CIRCULAR: True,
        OrbitType.ELLIPTICAL: True,
        OrbitType.PARABOLIC: False,
        OrbitType.HYPERBOLIC: False
    }
    return closure_dict[orbit_type]