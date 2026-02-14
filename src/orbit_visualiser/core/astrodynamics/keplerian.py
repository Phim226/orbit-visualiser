import numpy as np
from math import pi
from numpy.typing import NDArray
from typing import Callable, Literal
from orbit_visualiser.core.astrodynamics.types import OrbitType

def state_pf_from_e_rp(e: float, rp: float, mu: float, nu: float, state: Literal["pos", "vel", "both"] = "both") -> list[NDArray[np.float64]]:
    """
    Takes the eccentricity, radius of periapsis, gravitational parameter and true anomaly and returns
    the perifocal position and velocity vectors of a satellite at the given true anomaly on an
    analytical Keplerian orbit defined by the eccentricity and radius of periapsis.

    Parameters
    ----------
    e : float
        Eccentricity
    rp : float
        Radius of periapsis (km)
    mu : float
        Gravitational parameter (km^3/s^2)
    nu : float
        True anomaly (rads)
    state : Literal["pos", "vel", "both"]
        String literal indicating the state vector(s) being returned, default = both
    Returns
    -------
    list[NDArray[np.float64]]
        A list containing the numpy arrays of the perifocal position (km) and perifocal velocity (km/s)
    """
    p = semi_parameter(e, rp)
    h = specific_ang_momentum(mu, p)

    r = perifocal_position_eq(e, p)(nu)
    v = perifocal_velocity_eq(e, mu, h)(nu)

    if state == "pos":
        return [r]
    elif state == "vel":
        return [v]
    elif state == "both":
        return [r, v]

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
    def _callable(nu: float, asymp_anomaly: float = np.nan) -> NDArray[np.float64]:
        pf_pos_eq: Callable = lambda nu: p*(1/(1 + e*np.cos(nu)))*np.array([np.cos(nu), np.sin(nu)])

        # If the true anomaly at the true anomaly of the asymptote then the satellite is at infinity,
        # but since this is returning the perifocal position then we need to put the appropriate sign
        # in front of the x and y infinities.
        if np.isclose(abs(nu), asymp_anomaly, atol = 0.0001, rtol = 0):
                # The true anomaly and the x and y values 'near' infinity are guaranteed to be
                # non-zero when the true anomaly is close to the true anomaly of the asymptote,
                # so we can safely use the sign function.
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

def semi_parameter(e: float, rp: float) -> float:
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

def semimajor_axis(orbit_type: OrbitType, e: float, rp: float) -> float:
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

def semiminor_axis(orbit_type: OrbitType, closed: bool, e: float, a: float) -> float:
    """
    Calculates the semi-minor axis using the eccentricity and semi-major axis.

    Parameters
    ----------
    orbit_type : OrbitType
        The orbit type enum
    closed: bool
        Boolean representing if the orbit is closed
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

    elif closed:
        return a*np.sqrt(1 - e**2)

    return np.inf

def apoapsis(orbit_type: OrbitType, e: float, a: float) -> float:
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

def asymptote_anomaly(closed: bool, e: float) -> float:
    """
    Calculate the asymptote of the true anomaly for open orbits using the eccentricity.

    Parameters
    ----------
    closed : bool
        Boolean representing if the orbit is closed
    e : float
        Eccentricity

    Returns
    -------
    float
        The true anomaly of the asymptote (rads)
    """
    if closed:
        return np.nan

    return np.arccos(-1/e)

def turning_angle(closed: bool, e: float) -> float:
    """
    Calculate the turning angle for open orbits using the eccentricity

    Parameters
    ----------
    closed : bool
        Boolean representing if the orbit is closed
    e : float
        Eccentricity

    Returns
    -------
    float
        The turning angle (rads)
    """
    if closed:
        return np.nan

    return 2*np.arcsin(1/e)

def aiming_radius(orbit_type: OrbitType, b: float) -> float:
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
        nu: float,
        mu: float,
        h: float,
        e: float,
        asymp_anomaly: float
) -> NDArray[np.float64]:
    """
    Calculates the radial-azimuthal velocity vector using the gravitational parameter, angular
    momentum, eccentricity and true anomaly.

    Parameters
    ----------
    orbit_type : OrbitType
        The orbit type enum
    nu : float
        True anomaly (rads)
    mu : float
        Gravitational parameter (km^3/s^2)
    h : float
        Specific angular momentum (km^2/s)
    e : float
        Eccentricity
    asymp_anomaly : float
        The asymptote of the free anomaly (rads)

    Returns
    -------
    NDArray[np.float64]
        The radial-azimuthal velocity vector [v_radial, v_azimuthal] (km/s)
    """
    v_rad = 0.0 if orbit_type is OrbitType.CIRCULAR else (mu/h)*e*np.sin(nu)
    v_azim = (0.0 if np.isclose(abs(nu), asymp_anomaly, atol = 0.0001, rtol = 0)
              else (mu/h)*(1 + e*np.cos(nu)))

    return np.array([v_rad, v_azim])

def velocity_magnitude(velocity: NDArray[np.float64]) -> float:
    """
    Velocity magnitude helper function.

    Parameters
    ----------
    velocity : NDArray[np.float64]
        Velocity vector (km/s)

    Returns
    -------
    float
        The magnitude of the velocity (km/s)
    """
    return np.linalg.norm(velocity)

def radius(nu: float, asymp_anomaly: float, p: float, e: float) -> float:
    """
    Calculates the magnitude of the radius using the orbit equation and the semi-parameter.

    Parameters
    ----------
    nu : float
        True anomaly (rads)
    asymp_anomaly : float
        The asymptote of the free anomaly (rads)
    p : float
        Semi-parameter (km)
    e : float
        Eccentricity

    Returns
    -------
    float
        The length of the radius vector (km)
    """
    if np.isclose(abs(nu), asymp_anomaly, atol = 0.0001, rtol = 0):
        return np.inf

    return p/(1 + e*np.cos(nu))

def specific_energy(orbit_type: OrbitType, mu: float, a: float) -> float:
    """
    Calculates the specific energy using the gravitational parameter and semi-major axis.

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

def escape_velocity(nu: float, asymp_anomaly: float, mu: float, r: float) -> float:
    """
    Calculates escape velocity at the given orbital radius and gravitational parameter.

    Parameters
    ----------
    nu : float
        The true anomaly (rads)
    asymp_anomaly : float
        The true anomaly of the asymptot (rads)
    mu : float
        Gravitational parameter (km^3/s^2)
    r : float
        Orbital radius (km)

    Returns
    -------
    float
        The escape velocity (km/s)
    """
    if np.isclose(abs(nu), asymp_anomaly, atol = 0.0001, rtol = 0):
        return 0.0

    return np.sqrt(2*mu/r)

def excess_velocity(closed: bool, mu: float, a: float) -> float:
    """
    Calculates the hyperbolic excess velocity (the velocity magnitude at infinity) for open orbits
    from the semimajor axis and gravitational parameter.

    Parameters
    ----------
    closed : bool
        Boolean representing if the orbit is closed
    mu : float
        Gravitational parameter (km^3/s^2)
    a : float
        Semi-major axis (km)

    Returns
    -------
    float
        The hyperbolic excess velocity (km/s)
    """
    if closed:
        return np.nan

    return np.sqrt(mu/a)

def flight_angle(nu: float, asymp_anomaly: float, e: float) -> float:
    """
    Calculates the flight angle of a satellite at the given true anomaly and eccentricity.

    Parameters
    ----------
    nu : float
        True anomaly (rads)
    asymp_anomaly : float
        The true anomaly of the asymptot (rads)
    e : float
        Eccentricity

    Returns
    -------
    float
        The satellite flight angle (rads)
    """
    if np.isclose(abs(nu), asymp_anomaly, atol = 0.0001, rtol = 0):
        return (nu/abs(nu))*pi/2

    return np.arctan2(e*np.sin(nu), 1 + e*np.cos(nu))

def characteristic_energy(orbit_type: OrbitType, closed: bool, mu: float, a: float) -> float:
    """
    Calculates the characteristic energy from the semi-major axis and gravitational parameter.

    Parameters
    ----------
    orbit_type : OrbitType
        The orbit type enum
    closed : bool
        Boolean representing if the orbit is closed
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
    if closed:
        c3 *= -1

    return c3

def orbital_period(closed: bool, mu: float, a: float) -> float:
    """
    Calculates the orbital period of closed orbits.

    Parameters
    ----------
    closed : bool
        Boolean representing if the orbit is closed
    mu : float
        Gravitational parameter (km^3/s^2)
    a : float
        Semi-major axis (km)

    Returns
    -------
    float
        The orbital period (s)
    """
    if not closed:
        return np.nan

    return (2*pi/np.sqrt(mu))*np.sqrt(a)**3

def mean_motion(orbit_type: OrbitType, period: float, mu: float, p: float, a: float) -> float:
    """
    Calculates the mean motion, using different formulae based on orbit type.

    Parameters
    ----------
    orbit_type : OrbitType
        The orbit type enum
    period : float
        The orbital period (s)
    mu : float
        Gravitational parameter (km^3/s^2)
    p : float
        Semi-parameter (km)
    a : float
        Semi-major axis (km)

    Returns
    -------
    float
        The mean motion (rads/s)
    """
    if orbit_type in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL) :
        return 2*pi/period

    elif orbit_type is OrbitType.PARABOLIC:
        return 2*np.sqrt(mu/p)

    elif orbit_type is OrbitType.HYPERBOLIC:
        return 2*np.sqrt(mu/abs(a**3))

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


def time_since_periapsis(orbit_type: OrbitType, m_anomaly: float, period: float, p: float, h: float, e: float) -> float:
    """
    Calculates the time since periapsis, using different formulae based on orbit type.

    Parameters
    ----------
    orbit_type : OrbitType
        The orbit type enum
    m_anomaly : float
        The mean anomaly (rads)
    period : float
        The orbital period (s)
    p : float
        Semi-parameter (km)
    h : float
        Specific angular momentum (km^2/s)
    e : float
        Eccentricity

    Returns
    -------
    float
        Time since periapsis (s)
    """
    if np.isneginf(m_anomaly):
        return -np.inf

    elif np.isinf(m_anomaly):
        return np.inf

    if orbit_type in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
        return period*m_anomaly/(2*pi)

    elif orbit_type is OrbitType.PARABOLIC:
        return (p**2/h)*m_anomaly

    elif orbit_type is OrbitType.HYPERBOLIC:
        return (p**2/h)*(e**2 - 1)**(-1.5)*m_anomaly