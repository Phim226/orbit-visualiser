import pytest
import numpy as np
from numpy.typing import NDArray
from math import pi
from typing import Callable
from orbit_visualiser.core import Satellite, Orbit, CentralBody, asymptote_anomaly

# ---------- True anomaly grids --------------------
@pytest.fixture(scope = "session")
def closed_anomaly_grid() -> NDArray[np.float64]:
    return np.linspace(0, 2*pi, 10)

@pytest.fixture(scope = "session")
def open_anomaly_grid() -> Callable[[Orbit, int], NDArray[np.float64]]:
    def _create_grid(e: float, num: int = 20) -> NDArray[np.float64]:
        t_asymp = asymptote_anomaly(e)
        offset = 0.01
        return np.linspace(-t_asymp + offset, t_asymp - offset, num)
    return _create_grid

# ---------- Orbital object factories --------------
@pytest.fixture
def orbit_factory() -> Callable[[float, float], Orbit]:
    def _create(e: float = 0.0, rp: float = 50_000.0) -> Orbit:
        return Orbit(e, rp)
    return _create

@pytest.fixture
def central_body_factory() -> Callable[[float], CentralBody]:
    def _create(mu: float = 398_600.0) -> CentralBody:
        return CentralBody(mu)
    return _create

@pytest.fixture
def satellite_factory(
    orbit_factory: Callable[[float, float], Orbit],
    central_body_factory: Callable[[float], CentralBody]
) -> Callable[[float, float, float], Satellite]:
    def _create(e: float = 0.0, rp: float = 50_000.0, mu: float = 398_600.0) -> Satellite:
        orbit: Orbit = orbit_factory(e, rp)
        central_body: CentralBody = central_body_factory(mu)
        return Satellite(orbit, central_body)
    return _create