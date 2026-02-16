import itertools
from typing import Any

# Tests to implement
# TODO: Magnitude of angular momentum vector
# TODO: Magnitude of eccentricity vector
# TODO: Periapsis/apoapsis symmetry (check apoapsis has true anomaly pi ahead of periapsis in a closed orbit)
# TODO: Orbital period consistency ((contourpy) need to figure out how to parametrize over semi-major axis values and from that get appropriate values for e and rp to pass in the test objects)

def _create_test_cases(values_dict: dict[str, list[float]], tag_orbits: bool = True) -> list[list[Any]]:
    test_cases = [list(value) for value in itertools.product(*values_dict.values())]

    if "e" in values_dict and tag_orbits:
        for values in test_cases:
            e = values[0]

            if np.isclose(e, 0):
                orbit_type = OrbitType.CIRCULAR
            elif 0 < e < 1:
                orbit_type = OrbitType.ELLIPTICAL
            elif np.isclose(e, 1):
                orbit_type = OrbitType.PARABOLIC
            elif e > 1:
                orbit_type = OrbitType.HYPERBOLIC

            values.append(orbit_type)

    return test_cases

# Eccentricity test values for 0 <= e < 1
e_closed_typical_cases: list[float] = [0, 0.5]
e_closed_edge_cases: list[float] = [0.0001, 0.99999]
e_closed_test_cases: list[float] = e_closed_typical_cases + e_closed_edge_cases
e_closed_typical_copy = e_closed_typical_cases.copy()
e_closed_typical_copy.pop(0)
e_elliptical_test_cases: list[float] = e_closed_typical_copy + e_closed_edge_cases

# Eccentricity test values for e >= 1
e_open_typical_cases: list[float] = [1, 1.5]
e_open_edge_cases: list[float] = [1.0000001, 1000]
e_open_test_cases: list[float] = e_open_typical_cases + e_open_edge_cases
e_open_typical_copy = e_open_typical_cases.copy()
e_open_typical_copy.pop(0)
e_hyperbolic_test_cases: list[float] = e_open_typical_copy + e_open_edge_cases

# Aggregated eccentricity test values
e_typical_cases: list[float] = e_closed_typical_cases + e_open_typical_cases
e_edge_cases: list[float] = e_closed_edge_cases + e_open_edge_cases
e_test_cases: list[float] = e_closed_test_cases + e_open_test_cases

# Radius of periapsis test values
rp_typical_cases: list[float] = [6789, 20_000, 100_000, 1_000_000]
rp_edge_cases: list[float] = [10, 938_382_001_928_942_153]
rp_test_cases: list[float] = rp_typical_cases + rp_edge_cases

# Gravitation parameter test values
mu_typical_cases: list[float] = [4902.8, 398_600.0, 37_931_200.0]
mu_edge_cases: list[float] = [1.0]
mu_test_cases = mu_typical_cases + mu_edge_cases

# Test cases for variable eccentricity, radius of periapsis and gravitational parameter
full_test_cases = _create_test_cases({"e": e_test_cases, "rp": rp_test_cases, "mu": mu_test_cases})

# Test cases for variable eccentricity < 1, radius of periapsis and gravitational parameter
standard_closed_test_cases = _create_test_cases({"e": e_closed_test_cases, "rp": rp_test_cases, "mu": mu_test_cases})

# Test cases for variable eccentricity > 1, radius of periapsis and gravitational parameter
standard_open_test_cases = _create_test_cases({"e": e_open_test_cases, "rp": rp_test_cases, "mu": mu_test_cases})

# Test cases for constant eccentricity and variable radius of periapsis and gravitational parameter
rp_mu_test_cases = _create_test_cases({"rp": rp_test_cases, "mu": mu_test_cases})

# Test cases variable eccentricity with orbit type tags
e_tagged_test_cases = _create_test_cases({"e": e_test_cases})

# Test cases for typical values
typical_test_cases = _create_test_cases({"e": e_typical_cases, "rp": rp_typical_cases, "mu": mu_typical_cases})
typical_closed_test_cases = _create_test_cases({"e": e_closed_typical_cases, "rp": rp_typical_cases, "mu": mu_typical_cases}, False)



