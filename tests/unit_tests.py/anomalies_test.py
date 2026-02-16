import pytest
import numpy as np
from math import pi
from numpy.typing import NDArray
from orbit_visualiser.core import eccentric_anomaly, mean_anomaly, OrbitType
from tests.test_cases import e_elliptical_test_cases, e_hyperbolic_test_cases, e_closed_test_cases


def test_circular_eccentric_anomaly(closed_anomaly_grid: NDArray[np.float64]):
    """
    Tests that the eccentric anomaly is equivalent to the true anomaly for circular orbits.
    """
    nu = closed_anomaly_grid
    e_anomaly = eccentric_anomaly(OrbitType.CIRCULAR, 0, nu)
    assert np.allclose(nu, e_anomaly)

@pytest.mark.parametrize("e", e_elliptical_test_cases)
def test_elliptical_eccentric_anomaly(closed_anomaly_grid: NDArray[np.float64], e: float):
    """
    Tests that the eccentric anomaly for elliptical orbits is equivalent to the expected value.
    """
    nu = closed_anomaly_grid
    e_anomaly = eccentric_anomaly(OrbitType.ELLIPTICAL, e, nu)

    expected = 2*np.arctan(np.sqrt((1 - e)/(1 + e))*np.tan(nu/2))
    expected = np.where(expected < 0, expected + 2*pi, expected)

    assert np.allclose(e_anomaly, expected)

def test_parabolic_eccentric_anomaly(open_anomaly_grid: NDArray[np.float64]):
    """
    Tests that the eccentric anomaly for parabolic orbits is NaN.
    """
    nu = open_anomaly_grid(e = 1)
    e_anomaly = eccentric_anomaly(OrbitType.PARABOLIC, 1, nu)

    assert np.allclose(e_anomaly, np.full(len(nu), np.nan), equal_nan = True)

@pytest.mark.parametrize("e", e_hyperbolic_test_cases)
def test_hyperbolic_eccentric_anomaly(open_anomaly_grid: NDArray[np.float64], e: float):
    """
    Tests that the eccentric anomaly for hyperbolic orbits is equivalent to the expected value.
    """
    nu = open_anomaly_grid(e = e)
    e_anomaly = eccentric_anomaly(OrbitType.HYPERBOLIC, e, nu)

    expected = 2*np.arctanh(np.sqrt((e - 1)/(e + 1))*np.tan(nu/2))

    assert np.allclose(e_anomaly, expected)

@pytest.mark.parametrize("e", e_closed_test_cases)
def test_closed_mean_anomaly(closed_anomaly_grid: NDArray[np.float64], e: float):
    """
    Tests that the mean anomaly for closed orbits is equivalent to the expected value.
    """
    nu = closed_anomaly_grid
    e_anomaly = eccentric_anomaly(OrbitType.CIRCULAR, e, nu)
    m_anomaly = mean_anomaly(OrbitType.CIRCULAR, e, nu, e_anomaly)
    expected = e_anomaly - e*np.sin(e_anomaly)

    assert np.allclose(m_anomaly, expected)

def test_parabolic_mean_anomaly(closed_anomaly_grid: NDArray[np.float64]):
    """
    Tests that the mean anomaly for parabolic orbits is equivalent to the expected value.
    """
    nu = closed_anomaly_grid
    m_anomaly = mean_anomaly(OrbitType.PARABOLIC, 1, nu, np.nan)

    expected = 0.5*np.tan(nu/2) + (1/6)*(np.tan(nu/2))**3

    assert np.allclose(m_anomaly, expected)

@pytest.mark.parametrize("e", e_hyperbolic_test_cases)
def test_hyperbolic_mean_anomaly(open_anomaly_grid: NDArray[np.float64], e: float):
    """
    Tests that the mean anomaly for hyperbolic orbits is equivalent to the expected value.
    """
    nu = open_anomaly_grid(e = e)
    e_anomaly = eccentric_anomaly(OrbitType.HYPERBOLIC, e, nu)
    m_anomaly = mean_anomaly(OrbitType.HYPERBOLIC, e, nu, e_anomaly)
    expected = e*np.sinh(e_anomaly) - e_anomaly

    assert np.allclose(m_anomaly, expected)