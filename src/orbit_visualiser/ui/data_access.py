from typing import Literal
import numpy as np
from numpy.typing import NDArray
from math import pi
from orbit_visualiser.core import Satellite, OrbitType, perifocal_position, perifocal_to_eci_trans_mat

# TODO: Write docstrings
class OrbitDataAccess():

    def __init__(self, satellite: Satellite):
        self._satellite = satellite

    @property
    def satellite(self) -> Satellite:
        return self._satellite

    def _get_anomaly_data(self, num_points: int) -> NDArray[np.float64]:
        if self._satellite.orbit.orbit_type in (OrbitType.CIRCULAR, OrbitType.ELLIPTICAL):
            anomaly_range = (0, 2*pi)

        else:
            # If the orbit is open then we need a small offset so the plot doesn't evaluate to infinity
            # and cause runtime errors or unusual graphical artifacts.
            offset = 0.0001
            asymptote_anomaly = self._satellite.orbit.asymptote_anomaly
            anomaly_range = (-(asymptote_anomaly - offset), asymptote_anomaly - offset)

        return np.linspace(anomaly_range[0], anomaly_range[1], num_points)

    def get_orbit_data(self, num_points: int, ref_frame: Literal["eci", "perifocal"] = "eci") -> NDArray[np.float64]:
        orb = self._satellite.orbit
        pf_pos_data = perifocal_position(orb.eccentricity, orb.semi_parameter, self._get_anomaly_data(num_points))

        if ref_frame == "perifocal":
            return pf_pos_data

        eci_trans = perifocal_to_eci_trans_mat(orb.right_ascen_of_ascend_node, orb.inclination,
                                                    orb.argument_of_periapsis)
        return np.matmul(eci_trans, pf_pos_data)

    def get_sat_position(self, ref_frame: Literal["eci", "perifocal"] = "eci") -> NDArray[np.float64]:
        if ref_frame == "perifocal":
            pass

        return self._satellite.position