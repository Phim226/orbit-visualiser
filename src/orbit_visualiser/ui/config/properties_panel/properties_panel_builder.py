from tkinter import Tk
import numpy as np
from orbit_visualiser.ui import OrbitFigure, PropertySpec, VariableSpec
from orbit_visualiser.core import Orbit, Satellite, CentralBody


class PropertiesBuilder():


    def __init__(
            self, root: Tk,
            config_frame_placement: tuple[str],
            orbit_fig: OrbitFigure,
            orbit: Orbit,
            central_body: CentralBody,
            satellite: Satellite
    ):
        self._root = root

        self._orbit_fig = orbit_fig
        self._orbit = orbit
        self._central_body = central_body
        self._sat = satellite

        self._orbital_properties: dict[str, PropertySpec[Orbit]] = {
            "orbit_type" : PropertySpec("Orbit type", orbit, None, lambda orbit: orbit.orbit_type),
            "semi_major_axis" : PropertySpec("Semi-major axis", orbit, "km", lambda orbit: orbit.a),
            "semi_minor_axis" : PropertySpec("Semi-minor axis", orbit, "km", lambda orbit: orbit.b),
            "radius_apoapsis": PropertySpec("Radius of apoapsis", orbit, "km", lambda orbit: orbit.ra),
            "semi_parameter" : PropertySpec("Semi-parameter", orbit, "km", lambda orbit: orbit.p),
            "asymptote_anomal" : PropertySpec("Anomaly of asymptote", orbit, "°", lambda orbit: np.degrees(orbit.t_asymp)),
            "turn_angle" : PropertySpec("Turning angle", orbit, "°", lambda orbit: np.degrees(orbit.turn_angle)),
            "aim_radius" : PropertySpec("Aiming radius", orbit, "km", lambda orbit: orbit.aim_rad)
        }

        self._satellite_properties: dict[str, PropertySpec[Satellite]] = {
            "radius" : PropertySpec("Radius", satellite, "km", lambda sat: sat.r),
            "x_pos" : PropertySpec("Perifocal x position", satellite, "km", lambda sat: sat.pos_pf[0]),
            "y_pos" : PropertySpec("Perifocal y position", satellite, "km", lambda sat: sat.pos_pf[1]),
            "period" : PropertySpec("Orbital period", satellite, "s", lambda sat: sat.period),
            "mean_motion" : PropertySpec("Mean motion", satellite, "°/s", lambda sat: np.degrees(sat.n)),
            "e_anomaly" : PropertySpec("Eccentric anomaly", satellite, "°", lambda sat: np.degrees(sat.e_anomaly)),
            "m_anomaly" : PropertySpec("Mean anomaly", satellite, "°", lambda sat: np.degrees(sat.m_anomaly)),
            "time_periapsis" : PropertySpec("Time since periapsis", satellite, "s", lambda sat: sat.t_p),
            "ang_momentum" : PropertySpec("Angular momentum", satellite, "km²/s", lambda sat: sat.h),
            "velocity" : PropertySpec("Velocity", satellite, "km/s", lambda sat: sat.v),
            "azim_velocity" : PropertySpec("Azimuthal velocity", satellite, "km/s", lambda sat: sat.v_azim),
            "radial_velocity" : PropertySpec("Radial velocity", satellite, "km/s", lambda sat: sat.v_radial),
            "esc_velocity" : PropertySpec("Escape velocity", satellite, "km/s", lambda sat: sat.v_esc),
            "excess_velocity" : PropertySpec("Excess velocity", satellite, "km/s", lambda sat: sat.v_inf),
            "flight_angle" : PropertySpec("Flight angle", satellite, "°", lambda sat: np.degrees(sat.gam)),
            "spec_energy" : PropertySpec("Specific energy", satellite, "km²/s²", lambda sat: sat.eps),
            "char_energy" : PropertySpec("Characteristic energy", satellite, "km²/s²", lambda sat: sat.c3)
        }

        self._property_specs: dict[str, PropertySpec] = self._orbital_properties | self._satellite_properties