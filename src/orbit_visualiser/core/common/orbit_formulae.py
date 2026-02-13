import numpy as np
from numpy.typing import NDArray
from typing import Callable
from orbit_visualiser.core.common.types import OrbitType


def perifocal_position_eq(e: float, p: float) -> Callable[[float], NDArray[np.float64]]:
    """
    Perifocal orbit equation factory.

    Parameters
    ----------
    e : float
        Eccentricity
    p : float
        Semi-parameter (km)

    Returns
    -------
    Callable[[float], NDArray[np.float64]]
        The perifocal orbit equation, taking the true anomaly as an argument
    """
    def _callable(nu: float) -> NDArray[np.float64]:
        return p*(1/(1 + e*np.cos(nu)))*np.array([np.cos(nu), np.sin(nu)])
    return _callable

def perifocal_velocity_eq(e: float, mu: float, h: float) -> Callable[[float], NDArray[np.float64]]:
    """
    Perifocal velocity equation factory.

    Parameters
    ----------
    e  : float
        Eccentricity
    mu : float
        Gravitational parameter (km^3/s^2)
    h  : float
        Specific angular momentum (km^2/s)

    Returns
    -------
    Callable[[float], NDArray[np.float64]]
        The perifocal velocity equation, taking the true anomaly as an argument
    """
    def _callable(nu: float) -> NDArray[np.float64]:
        return (mu/h)*np.array([-np.sin(nu), e + np.cos(nu)])
    return _callable


def orbital_type(e: float) -> OrbitType:
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

def orbital_closure(orbit_type: OrbitType) -> bool:
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