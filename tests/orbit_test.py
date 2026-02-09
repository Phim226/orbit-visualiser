import pytest
import numpy as np
from typing import Callable
from orbit_visualiser.core import Orbit
from tests.test_cases import e_tagged_test_cases

# TODO: implement sanity tests (use different formulae for the same value and check they are equal)

@pytest.mark.parametrize("e, closure, orbit_type", e_tagged_test_cases)
def test_orbit_type_sanity(
    orbit_factory: Callable[[float, float], Orbit],
    e: float,
    closure: str,
    orbit_type: str):
    """
    If the eccentricity is 0<= e < 1 then the orbit is closed (True) and either circular (e = 0) or
    elliptical (0 < e < 1). If the eccentricity is e >= 1 then the orbit is not closed (False) and
    is either parabolic (e = 1) or hyperbolic (e > 1).
    """
    orbit: Orbit = orbit_factory(e = e)

    test_orbit_closure = "closed" if orbit.is_closed else "open"
    test_orbit_type = orbit.orbit_type

    assert test_orbit_closure == closure and test_orbit_type == orbit_type