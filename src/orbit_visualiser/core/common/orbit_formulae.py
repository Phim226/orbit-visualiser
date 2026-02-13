import numpy as np
from numpy.typing import NDArray
from typing import Callable

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