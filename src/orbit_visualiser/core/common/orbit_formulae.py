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
        The perifocal orbit equation, taking the true anomaly (rads) and the true anomaly of the
        asymptote (rads) as arguments
    """
    def _callable(nu: float, asymp_anomaly: float) -> NDArray[np.float64]:
        pf_pos_eq: Callable = lambda nu: p*(1/(1 + e*np.cos(nu)))*np.array([np.cos(nu), np.sin(nu)])

        # If the true anomaly at the true anomaly of the asymptote then the satellite is at infinity,
        # but since this is returning the perifocal position then we need to put the appropriate sign
        # in front of the x and y infinities.
        if np.isclose(abs(nu), asymp_anomaly, atol = 0.0001, rtol = 0):
                # The true anomaly and the x and y values 'near' infinity are guaranteed to be
                # non-zero in this part of the code, so we can safely use the sign function.
                nu_offset = np.sign(nu)*np.deg2rad(0.01)
                nu_close_to_inf = nu - nu_offset
                x_close_to_inf, y_close_to_inf = pf_pos_eq(nu_close_to_inf)

                return np.array([np.sign(x_close_to_inf)*np.inf, np.sign(y_close_to_inf)*np.inf])

        return pf_pos_eq(nu)

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

def semi_parameter_from_periapsis(e: float, rp: float) -> float:
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

def semimajor_axis_from_periapsis(orbit_type: OrbitType, e: float, rp: float) -> float:
    """
    Calculates the semimajor axis using the eccentricity and radius of periapsis.

    Parameters
    ----------
    orbit_type : OrbitType
        The orbit type enum
    e : float
        Eccentricity
    rp : float
        Radius of periapsis (km)

    Returns
    -------
    float
        Semi-major axis (km)
    """
    if orbit_type is OrbitType.PARABOLIC:
        return np.inf

    return rp/(1 - e)

def semiminor_axis_from_semimajor(orbit_type: OrbitType, e: float, a: float) -> float:
    """
    Calculates the semi-minor axis using the eccentricity and semi-major axis.

    Parameters
    ----------
    orbit_type : OrbitType
        The orbit type enum
    e : float
        Eccentricity
    a : float
        Semi-major axis (km)

    Returns
    -------
    float
        Semi-minor axis (km)
    """
    if orbit_type is OrbitType.HYPERBOLIC:
        return a*np.sqrt(e**2 - 1)

    elif orbit_closed(orbit_type):
        return a*np.sqrt(1 - e**2)

    return np.inf

def apoapsis_from_semimajor(orbit_type: OrbitType, e: float, a: float) -> float:
    """
    Calculates the radius of apoapsis using the eccentricity and semi-major axis.

    Parameters
    ----------
    orbit_type : OrbitType
        The orbit type enum
    e : float
        Eccentricity
    a : float
        Semi-major axis (km)

    Returns
    -------
    float
        Radius of apoapsis (km)
    """
    if orbit_type is OrbitType.PARABOLIC:
        return np.inf

    return a*(1 + e)

def asymptote_anomaly(orbit_type: OrbitType, e: float) -> float:
    """
    Calculate the asymptote of the true anomaly for open orbits using the eccentricity.

    Parameters
    ----------
    orbit_type : OrbitType
        The orbit type enum
    e : float
        Eccentricity

    Returns
    -------
    float
        The true anomaly of the asymptote (rads)
    """
    if orbit_closed(orbit_type):
        return np.nan

    return np.arccos(-1/e)

def turning_angle(orbit_type: OrbitType, e: float) -> float:
    """
    Calculate the turning angle for open orbits using the eccentricity

    Parameters
    ----------
    orbit_type : OrbitType
        The orbit type enum
    e : float
        Eccentricity

    Returns
    -------
    float
        The turning angle (rads)
    """
    if orbit_closed(orbit_type):
        return np.nan

    return 2*np.arcsin(1/e)

def aiming_radius_from_semiminor(orbit_type: OrbitType, b: float) -> float:
    """
    Calculate the aiming radius for hyperbolic orbits using the semi-minor axis

    Parameters
    ----------
    orbit_type : OrbitType
        The orbit type enum
    b : float
        Semi-minor axis (km)

    Returns
    -------
    float
        Aiming radius (km)
    """
    if orbit_type is not OrbitType.HYPERBOLIC:
        return np.nan

    return -b

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

def radial_azimuthal_velocity(
        orbit_type: OrbitType,
        nu: float
) -> NDArray[np.float64]:
    pass