from math import pi
from pytest import Subtests
import numpy as np
from orbit_visualiser.core import CentralBody, Orbit, Satellite


def test_circular_orbit_flight_angle(subtests: Subtests):
    satellite = Satellite(Orbit(), CentralBody())
    test_cases = np.linspace(0, 2*pi, num = 20)

    for i, nu in enumerate(test_cases):
        with subtests.test("Circular orbit flight angle test cases", i = i):
            satellite.nu = nu
            satellite.update_satellite_properties()

            gam = satellite.gam
            assert np.isclose(gam, 0)

def test_circular_orbit_radial_velocity(subtests: Subtests):
    satellite = Satellite(Orbit(), CentralBody())
    test_cases = np.linspace(0, 2*pi, num = 20)

    for i, nu in enumerate(test_cases):
        with subtests.test("Circular orbit radial velocity test cases", i = i):
            satellite.nu = nu
            satellite.update_satellite_properties()

            v_rad = satellite.v_radial
            assert np.isclose(v_rad, 0)

def test_circular_orbit_anomalies(subtests: Subtests):
    satellite = Satellite(Orbit(), CentralBody())
    test_cases = np.linspace(0, 2*pi, num = 20)

    for i, nu in enumerate(test_cases):
        with subtests.test("Circular orbit mean/eccentric/true anomaly test cases", i = i):
            satellite.nu = nu
            satellite.update_satellite_properties()

            mean = satellite.m_anomaly
            eccentric = satellite.e_anomaly
            assert np.isclose(mean, nu) and np.isclose(eccentric, nu)