import pytest
import numpy as np
from numpy.typing import NDArray
from typing import Callable
from orbit_visualiser.core import Orbit, Satellite
from tests.test_cases import rp_test_cases, rp_mu_test_cases

@pytest.mark.parametrize("rp", rp_test_cases)
def test_circular_orbit_distance_parameters(
    orbit_factory_from_elements: Callable[[float, float, float, float], Orbit],
    rp: float
):
    """
    For circular orbits the radius of periapsis, radius of apoapsis, semi-major axis, semi-minor axis
    and orbital parameter should all be equivalent.
    """
    orbit: Orbit = orbit_factory_from_elements(rp = rp)

    periapsis = orbit.radius_of_periapsis
    distances = [
        orbit.radius_of_apoapsis,
        orbit.semimajor_axis,
        orbit.semiminor_axis,
        orbit.semi_parameter
    ]
    assert np.allclose(distances, periapsis)

@pytest.mark.parametrize("rp, mu", rp_mu_test_cases)
def test_circular_orbit_flight_angle(
    satellite_factory_from_elements:  Callable[[float, float, float, float], Satellite],
    closed_anomaly_grid: NDArray[np.float64],
    rp: float,
    mu: float):
    """
    For circular orbits the flight angle of a satellite at any true anomaly should be 0.
    """


    for nu in closed_anomaly_grid:
        satellite: Satellite = satellite_factory_from_elements(rp = rp, mu = mu, nu = nu)

        gam = satellite.flight_angle
        assert np.isclose(gam, 0)

@pytest.mark.parametrize("rp, mu", rp_mu_test_cases)
def test_circular_orbit_radial_velocity(
    satellite_factory_from_elements:  Callable[[float, float, float, float], Satellite],
    closed_anomaly_grid: NDArray[np.float64],
    rp: float,
    mu: float):
    """
    For circular orbits the radial velocity of a satellite at any true anomaly should be 0.
    """


    for nu in closed_anomaly_grid:
        satellite: Satellite = satellite_factory_from_elements(rp = rp, mu = mu, nu = nu)

        v_rad = satellite.radial_azimuthal_velocity[0]
        assert np.isclose(v_rad, 0)

@pytest.mark.parametrize("rp, mu", rp_mu_test_cases)
def test_circular_orbit_anomalies(
    satellite_factory_from_elements:  Callable[[float, float, float, float], Satellite],
    closed_anomaly_grid: NDArray[np.float64],
    rp: float,
    mu: float):
    """
    For circular orbits the eccentric, mean and true anomalies should always be equivalent.
    """


    for nu in closed_anomaly_grid:
        satellite: Satellite = satellite_factory_from_elements(rp = rp, mu = mu, nu = nu)

        anomalies = [satellite.mean_anomaly, satellite.eccentric_anomaly]
        assert np.allclose(anomalies, nu)