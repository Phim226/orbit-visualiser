import pytest
from orbit_visualiser.core import OrbitType, orbit_type
from tests.test_cases import e_tagged_test_cases

@pytest.mark.parametrize("e, type", e_tagged_test_cases)
def test_classification(e: float, type: OrbitType):
    """
    Tests that the eccentricity is correctly classified as representing a circular, elliptical,
    parabolic or hyperbolic orbit.
    """
    assert orbit_type(e) is type