import itertools
from typing import Any

# Tests to implement
# TODO: Magnitude of angular momentum vector
# TODO: Magnitude of eccentricity vector
# TODO: Periapsis/apoapsis symmetry (check apoapsis has true anomaly pi ahead of periapsis in a closed orbit)
# TODO: Orbital period consistency (need to figure out how to parametrize over semi-major axis values and from that get appropriate values for e and rp to pass in the test objects)

def _create_test_cases(values_dict: dict[str, list[float]]) -> list[list[Any]]:
    test_cases = [list(value) for value in itertools.product(*values_dict.values())]

    if "e" in values_dict:
        for values in test_cases:
            e = values[0]
            values.append("closed" if e < 1 else "open")

            orbit_type: str
            if e == 0:
                orbit_type = "circular"
            elif 0 < e < 1:
                orbit_type = "elliptical"
            elif e == 1:
                orbit_type = "parabolic"
            else:
                orbit_type = "hyperbolic"
            values.append(orbit_type)


    return test_cases

e_test_cases: list[float] = [0, 0.0001, 0.1, 0.5, 0.9, 0.99999, 1, 1.0000001, 1.5, 2, 10, 1000]
rp_test_cases: list[float] = [10, 6789, 20_000, 50_000, 100_000, 1_000_000, 938_382_001_928_942_153]
mu_test_cases: list[float] = [1.0, 4902.8, 42_828.0, 398600.0, 37_931_200.0]

# Test cases for variable eccentricity, radius of periapsis and gravitational parameter
standard_test_cases = _create_test_cases({"e": e_test_cases, "rp": rp_test_cases, "mu": mu_test_cases})

# Test cases for constant eccentricity and variable radius of periapsis and gravitational parameter
rp_mu_test_cases = _create_test_cases({"rp": rp_test_cases, "mu": mu_test_cases})

# Test cases variable eccentricity with orbit type tags
e_tagged_test_cases = _create_test_cases({"e": e_test_cases})



