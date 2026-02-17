import pytest
from math import pi
import numpy as np
from numpy.typing import NDArray
from orbit_visualiser.core import (OrbitType, eccentricity_vector_from_state, eccentricity_from_state,
                                   true_anomaly_from_state, semi_parameter_from_momentum,
                                   semi_parameter_from_eccentricity, semimajor_axis,
                                   semiminor_axis, periapsis, apoapsis,
                                   asymptote_anomaly, turning_angle, aiming_radius,
                                   orbital_period, mean_motion)

@pytest.mark.parametrize("r, v, mu, expected", [
    (np.array([99_650.0, 0]), np.array([0, 2.0]), 398_600.0, np.array([0.0, 0.0])),
    (np.array([49_825.0, 0]), np.array([1.0, 2.0]), 398_600.0, np.array([-0.5, -0.25]))
])
def test_eccentricity_vector_from_state(
        r: NDArray[np.float64],
        v: NDArray[np.float64],
        mu: float,
        expected: NDArray[np.float64]
):
    """
    Test that the formula calculating the eccentricity vector from state gives the expected vector.
    """
    result = eccentricity_vector_from_state(r, v, mu)
    assert np.allclose(result, expected)

@pytest.mark.parametrize("r, v, mu, expected_mag", [
    (np.array([49_825.0, 0]), np.array([4.0, 2.0]), 398_600.0, 1.0)
])
def test_eccentricity_from_state(
        r: NDArray[np.float64],
        v: NDArray[np.float64],
        mu: float,
        expected_mag: float
):
    """
    Test that the formula calculating the eccentricity is in the correct magnitude range.
    """
    result = eccentricity_from_state(r, v, mu)
    assert result > expected_mag

@pytest.mark.parametrize("r, expected", [
    ([50_000.0, 0], 0),
    ([0, 50_000.0], pi/2),
    ([-50_000.0, 0], pi),
    ([0, -50_000.0], 3*pi/2)
])
def test_true_anomaly_from_state(r: NDArray[np.float64], expected: float):
    """
    Test that the formula for calculating the true anomaly from the current position gives the
    expected value.
    """
    result = true_anomaly_from_state(r)
    assert np.isclose(result, expected)

@pytest.mark.parametrize("h, mu, expected", [
    (100_000.0, 398600.0, 100_000**2/398600.0)
])
def test_semi_parameter_from_momentum(h: float, mu: float, expected: float):
    """
    Test that the formula for the semi-parameter using the specific angular momentum gives the
    expected value.
    """
    result = semi_parameter_from_momentum(h, mu)
    assert np.isclose(result, expected)

@pytest.mark.parametrize("e, rp, expected", [
    (0, 50_000, 50_000),
    (1, 50_000, 100_000)
])
def test_semi_parameter_from_eccentricity(e: float, rp: float, expected: float):
    """
    Test that the formula for the semi-parameter using the eccentricity and radius of periapsis
    gives the expected value.
    """
    result = semi_parameter_from_eccentricity(e, rp)
    assert np.isclose(result, expected)

@pytest.mark.parametrize("e, rp, expected", [
    (0, 50_000, 50_000),
    (1, 50_000, np.nan),
    (1.5, 50_000, -100_000)
])
def test_semimajor_axis(e: float, rp: float, expected: float):
    """
    Test that the formula for the semi-major axis using the eccentricity and radius of periapsis
    gives the expected value.
    """
    result = semimajor_axis(e, rp)
    assert np.isclose(result, expected, equal_nan = True)

@pytest.mark.parametrize("e, a, expected", [
    (0.6, 50_000, 40_000),
    (1, np.nan, np.nan),
    (13/12, -60_000, -25_000)
])
def test_semiminor_axis(e: float, a: float, expected: float):
    """
    Test that the formula for the semi-minor axis using the eccentricity and semi-major axis
    gives the expected value.
    """
    result = semiminor_axis(e, a)
    assert np.isclose(result, expected, equal_nan = True)

@pytest.mark.parametrize("p, e, expected", [
    (100_000, 1, 50_000)
])
def test_periapsis(p: float, e: float, expected: float):
    """
    Test that the formula for the periapsis using the eccentricity and semi-parameter gives the
    expected value.
    """
    result = periapsis(p, e)
    assert np.isclose(result, expected)

@pytest.mark.parametrize("e, a, expected", [
    (1, np.nan, np.nan),
    (0.5, 50_000, 75_000),
])
def test_apoapsis(e: float, a: float, expected: float):
    """
    Test that the formula for the periapsis using the eccentricity and semi-major axis gives the
    expected value.
    """
    result = apoapsis(e, a)
    assert np.isclose(result, expected, equal_nan = True)

@pytest.mark.parametrize("e, expected", [
    (0, np.nan),
    (0.5, np.nan),
    (1, pi)
])
def test_asymptote_anomaly(e: float, expected: float):
    """
    Test that the formula for the true anomaly of the asymptote using the eccentricity gives the expected
    value.
    """
    result = asymptote_anomaly(e)
    assert np.isclose(result, expected, equal_nan = True)

@pytest.mark.parametrize("e, expected", [
    (0, np.nan),
    (0.5, np.nan),
    (2, pi/3)
])
def test_turning_angle(e: float, expected: float):
    """
    Test that the formula for the turning angle using the eccentricity gives the expected value.
    """
    result = turning_angle(e)
    assert np.isclose(result, expected, equal_nan = True)

@pytest.mark.parametrize("b, expected", [
    (10_000, np.nan),
    (np.nan, np.nan),
    (-10_000, 10_000)
])
def test_aiming_radius(b: float, expected: float):
    """
    Test that the formula for the aiming radius using the semi-minor axis gives the expected value.
    """
    result = aiming_radius(b)
    assert np.isclose(result, expected, equal_nan = True)

@pytest.mark.parametrize("orbit_type, mu, a, expected", [
    (OrbitType.PARABOLIC, 398_600.0, 20_000, np.nan),
    (OrbitType.HYPERBOLIC, 398_600.0, 20_000, np.nan),
    (OrbitType.ELLIPTICAL, 398_600.0, 20_000, (2*pi/np.sqrt(398_600.0))*np.sqrt(20_000)**3)
])
def test_orbital_period(orbit_type: OrbitType, mu: float, a: float, expected: float):
    """
    Test that the formula for the orbital period using the gravitational parameter and semi-major
    axis gives the expected value.
    """
    result = orbital_period(orbit_type, mu, a)
    assert np.isclose(result, expected, equal_nan = True)

@pytest.mark.parametrize("orbit_type, mu, p, a, expected", [
    (OrbitType.PARABOLIC, 398_600.0, 10_000, 20_000, 2*np.sqrt(398_600.0/(10_000**3))),
    (OrbitType.HYPERBOLIC, 398_600.0, 10_000, -20_000, np.sqrt(398_600.0/abs(-20_000**3))),
    (OrbitType.ELLIPTICAL, 398_600.0, 10_000, 20_000, np.sqrt(398_600.0/abs(20_000**3))),
])
def test_mean_motion(orbit_type: OrbitType, mu: float, p: float, a: float, expected: float):
    """
    Test that the formula for the mean motion using the gravitational parameter and semi-major
    axis gives the expected value.
    """
    result = mean_motion(orbit_type, mu, p, a)
    assert np.isclose(result, expected)