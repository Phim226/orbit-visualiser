import pytest
import numpy as np
from typing import Callable
from numpy.typing import NDArray
from orbit_visualiser.core import Orbit, Satellite, OrbitType
from tests.test_cases import full_test_cases


@pytest.mark.parametrize("e, rp, mu, orbit_type", full_test_cases)
def test_specific_energy_conservation(
    satellite_factory: Callable[[float, float, float], Satellite],
    closed_anomaly_grid: NDArray[np.float64],
    open_anomaly_grid: Callable[[Orbit, int], NDArray[np.float64]],
    e: float,
    rp: float,
    mu: float,
    orbit_type: str
):
    satellite: Satellite = satellite_factory(e, rp, mu)

    if orbit_type in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
        anomaly_grid = closed_anomaly_grid
    else:
        anomaly_grid = open_anomaly_grid(satellite._orbit)

    for nu in anomaly_grid:
        satellite.nu = nu
        satellite.update_satellite_properties()

        specific_energy = satellite.eps

        specific_kin_energy = 0.5*satellite.v**2
        specific_pot_energy = -mu/satellite.r
        vis_viva_energy = specific_kin_energy + specific_pot_energy

        assert np.isclose(specific_energy, vis_viva_energy)

@pytest.mark.parametrize("e, rp, mu, orbit_type", full_test_cases)
def test_specific_and_characteristic_energy_relation(
    satellite_factory: Callable[[float, float, float], Satellite],
    closed_anomaly_grid: NDArray[np.float64],
    open_anomaly_grid: Callable[[Orbit, int], NDArray[np.float64]],
    e: float,
    rp: float,
    mu: float,
    orbit_type: str
):
    satellite: Satellite = satellite_factory(e, rp, mu)

    if orbit_type in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
        anomaly_grid = closed_anomaly_grid
    else:
        anomaly_grid = open_anomaly_grid(satellite._orbit)

    for nu in anomaly_grid:
        satellite.nu = nu
        satellite.update_satellite_properties()

        specific_energy = satellite.eps
        characteristic_energy = satellite.c3

        assert np.isclose(2*specific_energy, characteristic_energy)
