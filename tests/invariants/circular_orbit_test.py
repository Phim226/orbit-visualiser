import pytest
import numpy as np
from numpy.typing import NDArray
from typing import Callable
from orbit_visualiser.core import Orbit, Satellite
from tests.test_cases import rp_test_cases, rp_mu_test_cases

@pytest.mark.parametrize("rp", rp_test_cases)
def test_circular_orbit_distance_parameters(orbit_factory, rp: float):
    """
    For circular orbits the radius of periapsis, radius of apoapsis, semi-major axis, semi-minor axis
    and orbital parameter should all be equivalent.
    """
    orbit: Orbit = orbit_factory(rp = rp)

    periapsis = orbit.rp
    distances = [
        orbit.ra,
        orbit.a,
        orbit.b,
        orbit.p
    ]
    assert np.allclose(distances, periapsis)

@pytest.mark.parametrize("rp, mu", rp_mu_test_cases)
def test_circular_orbit_flight_angle(
    satellite_factory: Callable[[float, float, float], Satellite],
    closed_anomaly_grid: NDArray[np.float64],
    rp: float,
    mu: float):
    """
    For circular orbits the flight angle of a satellite at any true anomaly should be 0.
    """
    satellite: Satellite = satellite_factory(rp = rp, mu = mu)

    for nu in closed_anomaly_grid:
        satellite.nu = nu
        satellite.update_satellite_properties()

        gam = satellite.gam
        assert np.isclose(gam, 0)

@pytest.mark.parametrize("rp, mu", rp_mu_test_cases)
def test_circular_orbit_radial_velocity(
    satellite_factory: Callable[[float, float, float], Satellite],
    closed_anomaly_grid: NDArray[np.float64],
    rp: float,
    mu: float):
    """
    For circular orbits the radial velocity of a satellite at any true anomaly should be 0.
    """
    satellite: Satellite = satellite_factory(rp = rp, mu = mu)

    for nu in closed_anomaly_grid:
        satellite.nu = nu
        satellite.update_satellite_properties()

        v_rad = satellite.v_radial
        assert np.isclose(v_rad, 0)

@pytest.mark.parametrize("rp, mu", rp_mu_test_cases)
def test_circular_orbit_anomalies(
    satellite_factory: Callable[[float, float, float], Satellite],
    closed_anomaly_grid: NDArray[np.float64],
    rp: float,
    mu: float):
    """
    For circular orbits the eccentric, mean and true anomalies should always be equivalent.
    """
    satellite: Satellite = satellite_factory(rp = rp, mu = mu)

    for nu in closed_anomaly_grid:
        satellite.nu = nu
        satellite.update_satellite_properties()

        anomalies = [satellite.m_anomaly, satellite.e_anomaly]
        assert np.allclose(anomalies, nu)