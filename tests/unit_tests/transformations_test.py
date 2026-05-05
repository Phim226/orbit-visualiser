import pytest
import numpy as np
from numpy.typing import NDArray
from orbit_visualiser.core import eci_to_perifocal_trans_mat, perifocal_to_eci_trans_mat

# TODO: Write more test cases

@pytest.mark.parametrize("raan, i, omega, expected", [
    (0.0, 0.0, 0.0, [[1.0, 0.0, 0.0],
                     [0.0, 1.0, 0.0],
                     [0.0, 0.0, 1.0]])
])
def test_eci_to_perifocal_trans_mat(raan: float, i: float, omega: float, expected: NDArray[np.float64]):
    result = eci_to_perifocal_trans_mat(raan, i, omega)
    assert np.allclose(result, expected)

@pytest.mark.parametrize("raan, i, omega, expected", [
    (0.0, 0.0, 0.0, [[1.0, 0.0, 0.0],
                     [0.0, 1.0, 0.0],
                     [0.0, 0.0, 1.0]])
])
def test_perifocal_to_eci_trans_mat(raan: float, i: float, omega: float, expected: NDArray[np.float64]):
    result = perifocal_to_eci_trans_mat(raan, i, omega)
    assert np.allclose(result, expected)