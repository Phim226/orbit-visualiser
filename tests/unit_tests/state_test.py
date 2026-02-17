import pytest
from math import pi
import numpy as np
from numpy.typing import NDArray
from orbit_visualiser.core import (perifocal_position, perifocal_velocity, radial_azimuthal_velocity,
                                   speed, radius_from_state, radius_from_orbit_eq, escape_velocity,
                                   flight_angle, time_since_periapsis)

@pytest.mark.parametrize("e, p, nu, expected", [
    (0.5, 75_000.0, pi/2, np.array([0.0, 75_000.0]))
])
def test_perifocal_position(e: float, p: float, nu: float, expected: NDArray[np.float64]):
    """
    Test that the perifocal position equation gives the expected value.
    """
    result = perifocal_position(e, p, nu)
    assert np.allclose(result, expected)

@pytest.mark.parametrize("e, mu, h, nu, expected", [
    (0.5, 398_600.0, 199_300, pi/6, np.array([-1.0, 1 + np.sqrt(3)]))
])
def test_perifocal_velocity(e: float, mu: float, h:float, nu: float, expected: NDArray[np.float64]):
    """
    Test that the perifocal velocity equation gives the expected value.
    """
    result = perifocal_velocity(e, mu, h, nu)
    assert np.allclose(result, expected)

@pytest.mark.parametrize("nu, mu, h, e, asymp_anomaly, expected", [
    (0.0, 398_600.0, 141_173.6519, 0.0, np.nan, np.array([0.0, 398_600.0/141_173.6519])),
    (pi/6, 398_600.0, 199_300.0, 0.5, np.nan, np.array([0.5, 2 + 0.5*np.sqrt(3)])),
    (2*pi/3, 398_600.0, 244_519.9378, 2.0, 2*pi/3, np.array([(398600.0/244_519.9378)*np.sqrt(3), 0.0]))
])
def test_radial_azimuthal_velocity(
        nu: float,
        mu: float,
        h: float,
        e: float,
        asymp_anomaly: float,
        expected: NDArray[np.float64]
):
    """
    Test that the radial-azimuthal velocity formulae give the expected values.
    """
    result = radial_azimuthal_velocity(nu, mu, h, e, asymp_anomaly)
    assert np.allclose(result, expected)

@pytest.mark.parametrize("v, expected", [
    (np.array([1.0, 2.0]), np.sqrt(5))
])
def test_speed(v: NDArray[np.float64], expected: float):
    """
    Test that the speed is properly calculated from the velocity.
    """
    result = speed(v)
    assert np.isclose(result, expected)

@pytest.mark.parametrize("r, expected", [
    (np.array([1.0, 2.0]), np.sqrt(5))
])
def test_radius_from_state(r: NDArray[np.float64], expected: float):
    """
    Test that the radius calculated from state gives the expected value.
    """
    result = radius_from_state(r)
    assert np.isclose(result, expected)

@pytest.mark.parametrize("nu, asymp_anomaly, p, e, expected", [
    (0.0, np.nan, 75_000.0, 0.5, 50_000.0),
    (pi, pi, 100_000.0, 1, np.nan)
])
def test_radius_from_orbit_eq(nu: float, asymp_anomaly: float, p: float, e: float, expected: float):
    """
    Test that the radius calculated using the orbit equation gives the expected value.
    """
    result = radius_from_orbit_eq(nu, asymp_anomaly, p, e)
    assert np.isclose(result, expected, equal_nan = True)

@pytest.mark.parametrize("nu, asymp_anomaly, mu, r, expected", [
    (0.0, np.nan, 398_600.0, 50_000.0, np.sqrt(2*398_600.0/50_000.0)),
    (pi, pi, 398_600.0, np.nan, 0.0)
])
def test_escape_velocity(nu: float, asymp_anomaly: float, mu: float, r: float, expected: float):
    """
    Test that the escape velocity calculation gives the expected value.
    """
    result = escape_velocity(nu, asymp_anomaly, mu, r)
    assert np.isclose(result, expected, equal_nan = True)

@pytest.mark.parametrize("nu, asymp_anomaly, e, expected", [
    (0.0, np.nan, 0.5, 0.0),
    (pi, pi, 1, pi/2),
    (2*pi/3, np.nan, 0.5, pi/6)
])
def test_flight_angle(nu: float, asymp_anomaly: float, e: float, expected: float):
    """
    Test that the flight_angle calculation gives the expected value.
    """
    result = flight_angle(nu, asymp_anomaly, e)
    assert np.isclose(result, expected, equal_nan = True)

@pytest.mark.parametrize("m_anomaly, period, p, h, e, expected", [
    (pi/6, 111_266.9613, 50_000, 141_173.6519, 0.0, 111_266.9613*(pi/6)/(2*pi)),
    (pi/6, np.nan, 100_000, 199_649.6932, 1.0, ((100_000)**2/199_649.6932)*(pi/6)),
    (pi/6, np.nan, 125_000, 223_215.1429, 1.5, ((125_000)**2/223_215.1429)*(1.5**2 - 1)**(-1.5)*(pi/6))
])
def test_time_since_periapsis(
    m_anomaly: float,
    period: float,
    p: float,
    h: float,
    e: float,
    expected: float):
    """
    Test that the time_since_periapsis calculation gives the expected value.
    """
    result = time_since_periapsis(m_anomaly, period, p, h, e)
    assert np.isclose(result, expected, equal_nan = True)