import pytest
import numpy as np
from typing import Callable
from numpy.typing import NDArray
from orbit_visualiser.core import OrbitType, NewSatellite
from tests.test_cases import full_test_cases, standard_open_test_cases

# TODO: implement sanity tests (use different formulae for the same value and check they are equal)

@pytest.mark.parametrize("e, rp, mu, orbit_type", standard_open_test_cases)
def test_satellite_orbital_speed_sanity(
    satellite_factory_from_elements: Callable[[float, float, float, float], NewSatellite],
    open_anomaly_grid: Callable[[float, int], NDArray[np.float64]],
    e: float,
    rp: float,
    mu: float,
    orbit_type: str
):
    """
    Sanity test to check that different formulae for the orbital speed give the same result.
    """
    anomaly_grid = open_anomaly_grid(e = e)

    for nu in anomaly_grid:
        satellite: NewSatellite = satellite_factory_from_elements(e = e, rp = rp, mu = mu, nu = nu)

        v_vals = [
            np.hypot(satellite.orbit.hyperbolic_excess_velocity, satellite.escape_velocity),
            np.linalg.norm(satellite.radial_azimuthal_velocity)
        ]

        assert np.allclose(v_vals, satellite.speed)

@pytest.mark.parametrize("e, rp, mu, orbit_type", full_test_cases)
def test_satellite_azimuthal_velocity_sanity(
    satellite_factory_from_elements: Callable[[float, float, float, float], NewSatellite],
    closed_anomaly_grid: NDArray[np.float64],
    open_anomaly_grid: Callable[[float, int], NDArray[np.float64]],
    e: float,
    rp: float,
    mu: float,
    orbit_type: str
):
    """
    Sanity test to check that different formulae for the azimuthal velocity give the same result.
    """
    if orbit_type in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
        anomaly_grid = closed_anomaly_grid
    else:
        anomaly_grid = open_anomaly_grid(e)

    for nu in anomaly_grid:
        satellite: NewSatellite = satellite_factory_from_elements(e = e, rp = rp, mu = mu, nu = nu)

        v_azim_vals = [
            satellite.specific_angular_momentum/satellite.radius,
            satellite.speed*np.cos(satellite.flight_angle)
        ]


        assert np.allclose(v_azim_vals, satellite.radial_azimuthal_velocity[1])

@pytest.mark.parametrize("e, rp, mu, orbit_type", full_test_cases)
def test_satellite_radial_velocity_sanity(
    satellite_factory_from_elements: Callable[[float, float, float, float], NewSatellite],
    closed_anomaly_grid: NDArray[np.float64],
    open_anomaly_grid: Callable[[float, int], NDArray[np.float64]],
    e: float,
    rp: float,
    mu: float,
    orbit_type: str
):
    """
    Sanity test to check that different formulae for the radial velocity give the same result.
    """
    if orbit_type in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
        anomaly_grid = closed_anomaly_grid
    else:
        anomaly_grid = open_anomaly_grid(e)

    for nu in anomaly_grid:
        satellite: NewSatellite = satellite_factory_from_elements(e = e, rp = rp, mu = mu)

        v_radial = satellite.speed*np.sin(satellite.flight_angle)

        assert np.isclose(v_radial, satellite.radial_azimuthal_velocity[0])

