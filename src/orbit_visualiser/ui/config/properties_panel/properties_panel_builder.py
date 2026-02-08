from tkinter import Frame, Label, StringVar, LabelFrame
from typing import Callable
import numpy as np
from orbit_visualiser.core import Orbit, Satellite
from orbit_visualiser.ui.common.builder import Builder
from orbit_visualiser.ui.common.specs import PropertySpec


class PropertiesBuilder(Builder):


    def __init__(
            self,
            config_frame: Frame,
            orbit: Orbit,
            satellite: Satellite
    ):
        self._config_frame = config_frame

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

    @property
    def property_specs(self) -> dict[str, PropertySpec]:
        return self._property_specs

    def build_properties_frame(self, format_value: Callable) -> None:
        props_frame = Frame(self._config_frame, padx = 2)
        self._properties_frame = props_frame

        self._build_separator(props_frame, "Properties")
        orbital_props_frame = LabelFrame(
            props_frame, bd = 1, relief = "sunken", text = "Orbit", font = self._subtitle_font
        )
        self._populate_properties(orbital_props_frame, self._orbital_properties, format_value)
        orbital_props_frame.pack(side = "top", anchor = "nw", pady = (2, 0))

        sat_props_frame = LabelFrame(
            props_frame, bd = 1, relief = "sunken", text = "Satellite",
            font = self._subtitle_font, width = 244
        )
        self._populate_properties(sat_props_frame, self._satellite_properties, format_value)
        sat_props_frame.pack(side = "top", anchor = "nw", pady = (2, 0), fill = "x")

        props_frame.pack(side = "top", anchor = "n", pady = (2, 0))

    def _populate_properties(
            self,
            frame: LabelFrame,
            properties: dict[str, PropertySpec],
            format_value: Callable
    ) -> None:
        for i, (key, spec) in enumerate(properties.items()):
            self._build_property_row(frame, key, spec, i, format_value)

        # Setting the weight allows the grid manager to stretch labels in _build_display into available space.
        frame.grid_columnconfigure(0, weight = 0)
        frame.grid_columnconfigure(1, weight = 1)

    def _build_property_row(
            self,
            frame: LabelFrame,
            property: str,
            spec: PropertySpec,
            row: int,
            format_value: Callable
    ) -> None:
        init_value = spec.getter(spec.obj)

        var = StringVar(value = format_value(init_value, spec.units))
        self.__setattr__(f"{property}_str", var)

        name_label = Label(frame, text = spec.label + ":", anchor = "w", font = self._slider_font)
        name_label.grid(row = row, column = 0, sticky = "w", padx = (0, 6))

        value_label = Label(frame, textvariable = var, anchor = "e", width = 13, font = self._slider_font)
        value_label.grid(row = row, column = 1, sticky = "ew", padx = (0, 6))