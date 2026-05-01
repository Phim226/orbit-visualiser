import pytest
from math import pi
import numpy as np
from numpy.typing import NDArray
from orbit_visualiser.core import (OrbitType, eccentricity_vector_from_state, eccentricity_from_state,
                                   true_anomaly, semi_parameter_from_momentum,
                                   semi_parameter_from_eccentricity, semimajor_axis,
                                   semiminor_axis, radius_of_periapsis, radius_of_apoapsis,
                                   asymptote_anomaly, turning_angle, aiming_radius,
                                   orbital_period, mean_motion, inclination, node_line,
                                   right_ascen_of_ascending_node, argument_of_periapsis)

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

@pytest.mark.parametrize("r, e, v_r, expected", [
    (np.array([50_000.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0]), 0.0, 0.0),
    (np.array([0.0, 50_000.0, 0.0]), np.array([1.0, 0.0, 0.0]), 1.0, pi/2),
    (np.array([-50_000.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0]), 0.0, pi),
    (np.array([0.0, -50_000.0, 0.0]), np.array([1.0, 0.0, 0.0]), -1.0, 3*pi/2)
])
def test_true_anomaly_from_state(r: NDArray[np.float64], e: NDArray[np.float64], v_r: float, expected: float):
    """
    Test that the formula for calculating the true anomaly from the current position, eccentricity
    and radial speed gives the expected value.
    """
    result = true_anomaly(r, e, v_r)
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
    result = radius_of_periapsis(p, e)
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
    result = radius_of_apoapsis(e, a)
    assert np.isclose(result, expected, equal_nan = True)

@pytest.mark.parametrize("e, expected", [
    (0, np.nan),
    (0.5, np.nan),
    (1, pi)
])
def test_asymptote_anomaly(e: float, expected: float):
    """
    Test that the formula for the true anomaly of the asymptote using the eccentricity gives the
    expected value.
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

@pytest.mark.parametrize("e, mu, a, expected", [
    (1.0, 398_600.0, 20_000, np.nan),
    (1.5, 398_600.0, 20_000, np.nan),
    (0.5, 398_600.0, 20_000, (2*pi/np.sqrt(398_600.0))*np.sqrt(20_000)**3)
])
def test_orbital_period(e: float, mu: float, a: float, expected: float):
    """
    Test that the formula for the orbital period using the gravitational parameter and semi-major
    axis gives the expected value.
    """
    result = orbital_period(e, mu, a)
    assert np.isclose(result, expected, equal_nan = True)

@pytest.mark.parametrize("e, mu, p, a, expected", [
    (1.0, 398_600.0, 10_000, 20_000, 2*np.sqrt(398_600.0/(10_000**3))),
    (1.5, 398_600.0, 10_000, -20_000, np.sqrt(398_600.0/abs(-20_000**3))),
    (0.5, 398_600.0, 10_000, 20_000, np.sqrt(398_600.0/abs(20_000**3))),
])
def test_mean_motion(e: float, mu: float, p: float, a: float, expected: float):
    """
    Test that the formula for the mean motion using the gravitational parameter and semi-major
    axis gives the expected value.
    """
    result = mean_motion(e, mu, p, a)
    assert np.isclose(result, expected)

@pytest.mark.parametrize("h, expected", [
    (np.array([0.0, 0.0, 1.0]), 0.0),
    (np.array([1.0, 0.0, 0.0]), pi/2),
    (np.array([0.0, 0.0, -1.0]), pi)
])
def test_inclination(h: NDArray[np.float64], expected: float):
    """
    Test that the formula for the inclination using the specific angular momentum gives the expected
    value.
    """
    result = inclination(h)
    assert np.isclose(result, expected)

@pytest.mark.parametrize("h, expected", [
    (np.array([0.0, 0.0, 1.0]), np.array([1.0, 0.0, 0.0])),
    (np.array([1.0, 2.0, 0.0]), np.array([-2.0, 1.0, 0.0]))
])
def test_node_line(h: NDArray[np.float64], expected: NDArray[np.float64]):
    """
    Test that the formula for the node line using the specific angular momentum gives the expected
    value.
    """
    result = node_line(h)
    assert np.allclose(result, expected)

@pytest.mark.parametrize("node_line, expected", [
    (np.array([1.0, 0.0, 0.0]), 0.0),
    (np.array([0.0, 1.0, 0.0]), pi/2),
    (np.array([-1.0, 0.0, 0.0]), pi),
    (np.array([0.0, -1.0, 0.0]), 3*pi/2)
])
def test_right_ascension_of_ascending_node(node_line: NDArray[np.float64], expected: float):
    """
    Test that the formula for the right ascension of the ascending node using the node line gives
    the expected value.
    """
    result = right_ascen_of_ascending_node(node_line)
    assert np.isclose(result, expected, equal_nan = True)

@pytest.mark.parametrize("node_line, e, expected", [
    (np.array([1.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0]), 0.0),
    (np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), pi/2),
    (np.array([1.0, 0.0, 0.0]), np.array([-1.0, 0.0, 0.0]), pi),
    (np.array([1.0, 0.0, 0.0]), np.array([0.0, -1.0, -1.0]), 3*pi/2)
])
def test_argument_of_periapsis(node_line: NDArray[np.float64], e: NDArray[np.float64], expected: float):
    """
    Test that the formula for the argument of periapsis using the node line and the eccentricity
    gives the expected value.
    """
    result = argument_of_periapsis(node_line, e)
    assert np.isclose(result, expected, equal_nan = True)