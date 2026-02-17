import pytest
import numpy as np
from typing import Callable
from orbit_visualiser.core import Orbit, run_orbit_prop, get_init_conditions_from_elements
from tests.test_cases import typical_closed_test_cases

# TODO: Improve testing tolerances across difference orbit sizes and eccentricities.

@pytest.mark.parametrize("e, rp, mu", typical_closed_test_cases)
def test_1p_propagation_vel(
        orbit_factory_from_elements: Callable[[float, float, float, float], Orbit],
        e: float,
        rp: float,
        mu: float
):
    """
    Test that the velocity of a propagated satellite after a single period from periapsis for
    typical closed orbits is within an acceptable tolerance to the analytical solution.
    """
    orbit: Orbit = orbit_factory_from_elements(e, rp, mu, 0.0)
    init_conditions = get_init_conditions_from_elements(orbit)

    prop = run_orbit_prop(orbit, init_conditions, orbit.orbital_period)

    vel_at_rp = init_conditions[2:]
    prop_vel_at_rp = prop.y[2:, -1]

    assert np.allclose(vel_at_rp, prop_vel_at_rp, atol = 1e-6)


@pytest.mark.parametrize("e, rp, mu", typical_closed_test_cases)
def test_1p_propagation_pos(
        orbit_factory_from_elements: Callable[[float, float, float, float], Orbit],
        e: float,
        rp: float,
        mu: float
):
    """
    Test that the position of a propagated satellite after a single period from periapsis for
    typical closed orbits is within an acceptable tolerance to the analytical solution.
    """
    orbit: Orbit = orbit_factory_from_elements(e, rp, mu, 0.0)
    init_conditions = get_init_conditions_from_elements(orbit)

    prop = run_orbit_prop(orbit, init_conditions, orbit.orbital_period)

    pos_at_rp = init_conditions[:2]
    prop_pos_at_rp = prop.y[:2, -1]

    assert np.allclose(pos_at_rp, prop_pos_at_rp, atol = 0.01)


@pytest.mark.parametrize("e, rp, mu", typical_closed_test_cases)
def test_1p_phase_shift(
        orbit_factory_from_elements: Callable[[float, float, float, float], Orbit],
        e: float,
        rp: float,
        mu: float
):
    """
    Test that the phase shift of a propagated satellite after a single period from periapsis for
    typical closed orbits is within an acceptable tolerance to the analytical solution.
    """
    orbit: Orbit = orbit_factory_from_elements(e, rp, mu, 0.0)
    init_conditions = get_init_conditions_from_elements(orbit)

    prop = run_orbit_prop(orbit, init_conditions, orbit.orbital_period)

    prop_pos_at_rp = prop.y[:2, -1]
    phase_shift = np.arctan2(prop_pos_at_rp[1], prop_pos_at_rp[0])

    assert np.isclose(phase_shift, 0, atol=1e-8)