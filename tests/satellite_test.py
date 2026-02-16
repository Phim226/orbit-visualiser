import pytest
import numpy as np
from typing import Callable
from numpy.typing import NDArray
from orbit_visualiser.core import Satellite, Orbit, OrbitType
from tests.test_cases import full_test_cases, standard_open_test_cases

# TODO: implement sanity tests (use different formulae for the same value and check they are equal)

@pytest.mark.parametrize("e, rp, mu, orbit_type", standard_open_test_cases)
def test_satellite_velocity_sanity(
    satellite_factory: Callable[[float, float, float], Satellite],
    open_anomaly_grid: Callable[[Orbit, int], NDArray[np.float64]],
    e: float,
    rp: float,
    mu: float,
    orbit_type: str
):
    """
    Sanity test to check that different formulae for the velocity give the same result.
    """
    satellite: Satellite = satellite_factory(e, rp, mu)

    anomaly_grid = open_anomaly_grid(satellite._orbit)

    for nu in anomaly_grid:
        satellite.nu = nu
        satellite.update_satellite_properties()

        v_vals = [
            np.hypot(satellite.v_inf, satellite.v_esc),
            np.hypot(satellite.v_azim, satellite.v_radial)
        ]

        assert np.allclose(v_vals, satellite.v)

@pytest.mark.parametrize("e, rp, mu, orbit_type", full_test_cases)
def test_satellite_azimuthal_velocity_sanity(
    satellite_factory: Callable[[float, float, float], Satellite],
    closed_anomaly_grid: NDArray[np.float64],
    open_anomaly_grid: Callable[[Orbit, int], NDArray[np.float64]],
    e: float,
    rp: float,
    mu: float,
    orbit_type: str
):
    """
    Sanity test to check that different formulae for the azimuthal velocity give the same result.
    """
    satellite: Satellite = satellite_factory(e, rp, mu)

    if orbit_type in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
        anomaly_grid = closed_anomaly_grid
    else:
        anomaly_grid = open_anomaly_grid(satellite._orbit)

    for nu in anomaly_grid:
        satellite.nu = nu
        satellite.update_satellite_properties()

        v_azim_vals = [
            satellite.h/satellite.r,
            satellite.v*np.cos(satellite.gam)
        ]


        assert np.allclose(v_azim_vals, satellite.v_azim)

@pytest.mark.parametrize("e, rp, mu, orbit_type", full_test_cases)
def test_satellite_radial_velocity_sanity(
    satellite_factory: Callable[[float, float, float], Satellite],
    closed_anomaly_grid: NDArray[np.float64],
    open_anomaly_grid: Callable[[Orbit, int], NDArray[np.float64]],
    e: float,
    rp: float,
    mu: float,
    orbit_type: str
):
    """
    Sanity test to check that different formulae for the radial velocity give the same result.
    """
    satellite: Satellite = satellite_factory(e, rp, mu)

    if orbit_type in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
        anomaly_grid = closed_anomaly_grid
    else:
        anomaly_grid = open_anomaly_grid(satellite._orbit)

    for nu in anomaly_grid:
        satellite.nu = nu
        satellite.update_satellite_properties()

        v_radial = satellite.v*np.sin(satellite.gam)

        assert np.isclose(v_radial, satellite.v_radial)

