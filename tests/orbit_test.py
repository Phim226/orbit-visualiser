from pytest import Subtests
import numpy as np
from orbit_visualiser.core import Orbit



def test_circular_orbit_distance_parameters(subtests: Subtests):
    circular_orbit_test_cases = [
        0.000000001,
        1,
        6789,
        50000,
        1000000,
        4095384905805829034
    ]

    for i, periapsis in enumerate(circular_orbit_test_cases):
        with subtests.test("Circular orbit distance parameters test cases", i = i):
            orbit = Orbit(rp = periapsis)
            rp = orbit.rp
            ra = orbit.ra
            a = orbit.a
            b = orbit.b
            p = orbit.p
            assert (np.isclose(rp, periapsis) and np.isclose(rp, ra) and
                    np.isclose(rp, a) and np.isclose(rp, b) and np.isclose(rp, p))

def test_orbit_type(subtests: Subtests):
    eccentricity_test_cases = {
        0 : ("circular", True),
        0.0000000001: ("elliptical", True),
        0.5: ("elliptical", True),
        0.9999999999: ("elliptical", True),
        1: ("parabolic", False),
        1.5: ("hyperbolic", False),
        10000000: ("hyperbolic", False)
    }

    for i, (e, output) in enumerate(eccentricity_test_cases.items()):
        with subtests.test("Orbit type test cases", i = i):
            orbit = Orbit(e = e)
            assert orbit.orbit_type == output[0] and orbit.is_closed is output[1]