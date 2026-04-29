import pytest
from orbit_visualiser.core import OrbitType, OrbitMotion, orbit_type, orbit_motion_type
from tests.test_cases import orbit_type_test_cases, motion_type_test_cases

@pytest.mark.parametrize("e, type", orbit_type_test_cases)
def test_orbit_type(e: float, type: OrbitType):
    """
    Tests that the eccentricity is correctly classified as representing a circular, elliptical,
    parabolic or hyperbolic orbit.
    """
    assert orbit_type(e) is type

@pytest.mark.parametrize("i, type", motion_type_test_cases)
def test_orbit_motion_type(i: float, type: OrbitMotion):
    """
    Tests that the inclination is correctly classified as prograde or retrograde.
    """
    assert orbit_motion_type(i) is type