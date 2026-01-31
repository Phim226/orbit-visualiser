from math import pi
from pytest import Subtests
import numpy as np
from orbit_visualiser.core import CentralBody, Orbit, Satellite


def test_circular_orbit_flight_angle(subtests: Subtests):
    """
    For circular orbits the flight angle of a satellite at any true anomaly should be 0.

    test_cases is the interval [0, 2pi] radians split into 20, so each value of the true
    anomaly increases by roughly 0.314 radians for each test.
    """
    satellite = Satellite(Orbit(), CentralBody())
    test_cases = np.linspace(0, 2*pi, num = 20)

    for i, nu in enumerate(test_cases):
        with subtests.test("Circular orbit flight angle test cases", i = i):
            satellite.nu = nu
            satellite.update_satellite_properties()

            gam = satellite.gam
            assert np.isclose(gam, 0)

def test_circular_orbit_radial_velocity(subtests: Subtests):
    """
    For circular orbits the radial velocity of a satellite at any true anomaly should be 0.

    test_cases is the interval [0, 2pi] radians split into 20, so each value of the true
    anomaly increases by roughly 0.314 radians for each test.
    """
    satellite = Satellite(Orbit(), CentralBody())
    test_cases = np.linspace(0, 2*pi, num = 20)

    for i, nu in enumerate(test_cases):
        with subtests.test("Circular orbit radial velocity test cases", i = i):
            satellite.nu = nu
            satellite.update_satellite_properties()

            v_rad = satellite.v_radial
            assert np.isclose(v_rad, 0)

def test_circular_orbit_anomalies(subtests: Subtests):
    """
    For circular orbits the eccentric, mean and true anomalies should always be equivalent.

    test_cases is the interval [0, 2pi] radians split into 20, so each value of the true
    anomaly increases by roughly 0.314 radians for each test.
    """
    satellite = Satellite(Orbit(), CentralBody())
    test_cases = np.linspace(0, 2*pi, num = 20)

    for i, nu in enumerate(test_cases):
        with subtests.test("Circular orbit mean/eccentric/true anomaly test cases", i = i):
            satellite.nu = nu
            satellite.update_satellite_properties()

            anomalies = [satellite.m_anomaly, satellite.e_anomaly]
            assert np.allclose(anomalies, nu)