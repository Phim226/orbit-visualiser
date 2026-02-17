import pytest
import numpy as np
from typing import Callable
from numpy.typing import NDArray
from orbit_visualiser.core import OrbitType, NewSatellite
from tests.test_cases import full_test_cases


@pytest.mark.parametrize("e, rp, mu, orbit_type", full_test_cases)
def test_specific_energy_conservation(
    satellite_factory_from_elements: Callable[[float, float, float, float], NewSatellite],
    closed_anomaly_grid: NDArray[np.float64],
    open_anomaly_grid: Callable[[float, int], NDArray[np.float64]],
    e: float,
    rp: float,
    mu: float,
    orbit_type: str
):
    """
    Test that the specific orbital energy remains constant across a given orbit.
    """

    if orbit_type in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
        anomaly_grid = closed_anomaly_grid
    else:
        anomaly_grid = open_anomaly_grid(e)

    for nu in anomaly_grid:
        satellite: NewSatellite = satellite_factory_from_elements(e = e, rp = rp, mu = mu, nu = nu)

        specific_energy = satellite.orbit.specific_energy

        specific_kin_energy = 0.5*satellite.speed**2
        specific_pot_energy = -mu/satellite.radius
        vis_viva_energy = specific_kin_energy + specific_pot_energy

        assert np.isclose(specific_energy, vis_viva_energy)

@pytest.mark.parametrize("e, rp, mu, orbit_type", full_test_cases)
def test_specific_and_characteristic_energy_relation(
    satellite_factory_from_elements: Callable[[float, float, float, float], NewSatellite],
    closed_anomaly_grid: NDArray[np.float64],
    open_anomaly_grid: Callable[[float, int], NDArray[np.float64]],
    e: float,
    rp: float,
    mu: float,
    orbit_type: str
):
    """
    Test that the characteristic energy is double the specific energy at various points across an orbit.
    """

    if orbit_type in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
        anomaly_grid = closed_anomaly_grid
    else:
        anomaly_grid = open_anomaly_grid(e)

    for nu in anomaly_grid:
        satellite: NewSatellite = satellite_factory_from_elements(e = e, rp = rp, mu = mu)

        specific_energy = satellite.orbit.specific_energy
        characteristic_energy = satellite.orbit.characteristic_energy

        assert np.isclose(2*specific_energy, characteristic_energy)
