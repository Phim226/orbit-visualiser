import pytest
import numpy as np
from typing import Callable
from numpy.typing import NDArray
from orbit_visualiser.core import (OrbitType, specific_ang_momentum_from_state, specific_ang_momentum,
                                   specific_orbital_energy, characteristic_energy, excess_velocity,
                                   vis_viva_speed)
from tests.test_cases import full_test_cases, standard_open_test_cases

@pytest.mark.parametrize("r, v, expected", [
    (np.array([1.0, 0]), np.array([0, 1.0]), np.array(1.0)),
    (np.array([1.0, 2.0]), np.array([2.0, 4.0]), np.array(0.0)),
    (np.array([50_000, 0]), np.array([0, 4.0]), np.array(4*50_000))
])
def test_specific_angular_momentum_from_state(
    r: NDArray[np.float64],
    v: NDArray[np.float64],
    expected: float
):
    """
    Test that the formula calculating the specific angular momentum from the state gives the
    expected value.
    """
    result = specific_ang_momentum_from_state(r, v)
    assert np.allclose(result, expected)

@pytest.mark.parametrize("mu, p, expected", [
    (398_600.0, 50_000.0, np.sqrt(398_600.0*50_000.0))
])
def test_specific_angular_momentum(mu: float, p: float, expected: float):
    """
    Test that the formula calculating the specific angular momentum from the gravitational parameter
    and semi-parameter gives the expected value.
    """
    result = specific_ang_momentum(mu, p)
    assert np.isclose(result, expected)

@pytest.mark.parametrize("orbit_type, mu, a, expected", [
    (OrbitType.PARABOLIC, 398_600.0, 50_000.0, 0.0),
    (OrbitType.CIRCULAR, 398_600.0, 50_000.0, -398_600.0/(2*50_000)),
    (OrbitType.HYPERBOLIC, 398_600.0, -50_000.0, 398_600.0/(2*50_000))
])
def test_specific_orbital_energy(orbit_type: OrbitType, mu: float, a: float, expected: float):
    """
    Test that the formula calculating the specific orbital energy gives the expected value.
    """
    result = specific_orbital_energy(orbit_type, mu, a)
    assert np.isclose(result, expected)

@pytest.mark.parametrize("orbit_type, mu, a, expected", [
    (OrbitType.PARABOLIC, 398_600.0, 50_000.0, 0.0),
    (OrbitType.CIRCULAR, 398_600.0, 50_000.0, -398_600.0/50_000),
    (OrbitType.HYPERBOLIC, 398_600.0, -50_000.0, 398_600.0/50_000)
])
def test_characteristic_energy(orbit_type: OrbitType, mu: float, a: float, expected: float):
    """
    Test that the formula calculating the characteristic energy gives the expected value.
    """
    result = characteristic_energy(orbit_type, mu, a)
    assert np.isclose(result, expected)

@pytest.mark.parametrize("orbit_type, mu, a, expected", [
    (OrbitType.CIRCULAR, 398_600.0, 50_000.0, np.nan),
    (OrbitType.HYPERBOLIC, 398_600.0, -50_000.0, np.sqrt(398_600.0/abs(-50_000.0)))
])
def test_excess_velocity(orbit_type: OrbitType, mu: float, a: float, expected: float):
    """
    Test that the formula calculating the excess velocity gives the expected value.
    """
    result = excess_velocity(orbit_type, mu, a)
    assert np.isclose(result, expected, equal_nan = True)

@pytest.mark.parametrize("r, a, mu, expected", [
    (20_000.0, 15_000.0, 398_600.0, np.sqrt(398_600.0*((2/20_000.0)-(1/15_000))))
])
def test_vis_viva_speed(r: float, a: float, mu: float, expected: float):
    """
    Test that the formula calculating the speed using the vis-viva equation gives the expected
    value.
    """
    result = vis_viva_speed(r, a, mu)
    assert np.isclose(result, expected, equal_nan = True)