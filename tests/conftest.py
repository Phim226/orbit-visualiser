import pytest
import numpy as np
from numpy.typing import NDArray
from math import pi
from typing import Callable
from orbit_visualiser.core import CentralBody, NewSatellite, NewOrbit, asymptote_anomaly

# ---------- True anomaly grids --------------------
@pytest.fixture(scope = "session")
def closed_anomaly_grid() -> NDArray[np.float64]:
    return np.linspace(0, 2*pi, 10)

@pytest.fixture(scope = "session")
def open_anomaly_grid() -> Callable[[float, int], NDArray[np.float64]]:
    def _create_grid(e: float, num: int = 20) -> NDArray[np.float64]:
        t_asymp = asymptote_anomaly(e)
        offset = 0.01
        return np.linspace(-t_asymp + offset, t_asymp - offset, num)
    return _create_grid

# ---------- Orbital object factories --------------
@pytest.fixture
def orbit_factory_from_state() -> Callable[[NDArray[np.float64], NDArray[np.float64], float], NewOrbit]:
    def _create(
            r: NDArray[np.float64] = np.array([50_000.0, 0]),
            v: NDArray[np.float64] = np.array([0, 2.823473039]),
            mu: float = 398_600.0
    ) -> NewOrbit:
        return NewOrbit(r, v, mu)
    return _create

@pytest.fixture
def orbit_factory_from_elements() -> Callable[[float, float, float, float], NewOrbit]:
    def _create(e: float = 0.0, rp: float = 50_000.0, mu: float = 398_600.0, nu: float = 0.0) -> NewOrbit:
        return NewOrbit.from_orbital_elements(e, rp, mu, nu)
    return _create

@pytest.fixture
def central_body_factory() -> Callable[[float], CentralBody]:
    def _create(mu: float = 398_600.0) -> CentralBody:
        return CentralBody(mu)
    return _create

@pytest.fixture
def satellite_factory_from_state(
    central_body_factory: Callable[[float], CentralBody]
) -> Callable[[NDArray[np.float64], NDArray[np.float64], float], NewSatellite]:
    def _create(
            r: NDArray[np.float64] = np.array([50_000.0, 0]),
            v: NDArray[np.float64] = np.array([0, 2.823473039]),
            mu: float = 398_600.0
    ) -> NewSatellite:
        central_body: CentralBody = central_body_factory(mu)
        return NewSatellite(r, v, central_body)
    return _create

@pytest.fixture
def satellite_factory_from_elements(
    orbit_factory_from_elements: Callable[[float, float, float, float], NewOrbit],
    central_body_factory: Callable[[float], CentralBody]
) -> Callable[[float, float, float], NewSatellite]:
    def _create(e: float = 0.0, rp: float = 50_000.0, mu: float = 398_600.0, nu: float = 0.0) -> NewSatellite:
        orbit: NewOrbit = orbit_factory_from_elements(e, rp, mu, nu)
        r = orbit.position
        v = orbit.velocity
        central_body: CentralBody = central_body_factory(mu)
        return NewSatellite(r, v, central_body)
    return _create