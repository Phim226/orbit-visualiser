from tkinter import Frame, Label, StringVar, LabelFrame
from typing import Callable
import numpy as np
from orbit_visualiser.core import Orbit, Satellite
from orbit_visualiser.ui.common.builder import Builder
from orbit_visualiser.ui.common.specs import PropertySpec
from orbit_visualiser.core.satellite import NewSatellite


class PropertiesBuilder(Builder):


    def __init__(
            self,
            config_frame: Frame,
            satellite: NewSatellite
    ):
        self._config_frame = config_frame
        self._satellite = satellite

        self._orbital_properties: dict[str, PropertySpec[NewSatellite]] = {
            "orbit_type" : PropertySpec("Orbit type", None, lambda sat: sat.orbit.orbit_type),
            "semi_major_axis" : PropertySpec("Semi-major axis", "km", lambda sat: sat.orbit.semimajor_axis),
            "semi_minor_axis" : PropertySpec("Semi-minor axis", "km", lambda sat: sat.orbit.semiminor_axis),
            "radius_apoapsis": PropertySpec("Radius of apoapsis", "km", lambda sat: sat.orbit.radius_of_apoapsis),
            "semi_parameter" : PropertySpec("Semi-parameter", "km", lambda sat: sat.orbit.semi_parameter),
            "period" : PropertySpec("Orbital period", "s", lambda sat: sat.orbit.orbital_period),
            "mean_motion" : PropertySpec("Mean motion", "°/s", lambda sat: np.degrees(sat.orbit.mean_motion)),
            "asymptote_anomal" : PropertySpec("Anomaly of asymptote", "°", lambda sat: np.degrees(sat.orbit.asymptote_anomaly)),
            "turn_angle" : PropertySpec("Turning angle", "°", lambda sat: np.degrees(sat.orbit.turning_angle)),
            "aim_radius" : PropertySpec("Aiming radius", "km", lambda sat: sat.orbit.aiming_radius),
            "spec_energy" : PropertySpec("Specific orbital energy", "km²/s²", lambda sat: sat.orbit.specific_orbital_energy),
            "char_energy" : PropertySpec("Characteristic energy", "km²/s²", lambda sat: sat.orbit.characteristic_energy),
            "excess_velocity" : PropertySpec("Excess velocity", "km/s", lambda sat: sat.orbit.hyperbolic_excess_velocity)
        }

        self._satellite_properties: dict[str, PropertySpec[NewSatellite]] = {
            "radius" : PropertySpec("Radius", "km", lambda sat: sat.radius),
            "x_pos" : PropertySpec("Perifocal x position", "km", lambda sat: sat.position[0]),
            "y_pos" : PropertySpec("Perifocal y position", "km", lambda sat: sat.position[1]),
            "e_anomaly" : PropertySpec("Eccentric anomaly", "°", lambda sat: np.degrees(sat.eccentric_anomaly)),
            "m_anomaly" : PropertySpec("Mean anomaly", "°", lambda sat: np.degrees(sat.mean_anomaly)),
            "time_periapsis" : PropertySpec("Time since periapsis", "s", lambda sat: sat.time_since_periapsis),
            "ang_momentum" : PropertySpec("Angular momentum", "km²/s", lambda sat: sat.specific_angular_momentum),
            "speed" : PropertySpec("Orbital speed", "km/s", lambda sat: sat.speed),
            "azim_velocity" : PropertySpec("Azimuthal velocity", "km/s", lambda sat: sat.radial_azimuthal_velocity[1]),
            "radial_velocity" : PropertySpec("Radial velocity", "km/s", lambda sat: sat.radial_azimuthal_velocity[0]),
            "flight_angle" : PropertySpec("Flight angle", "°", lambda sat: np.degrees(sat.flight_angle)),
            "esc_velocity" : PropertySpec("Escape velocity", "km/s", lambda sat: sat.escape_velocity)
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
            props_frame, bd = 1, relief = "sunken", text = "Current satellite state",
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
        init_value = spec.getter(self._satellite)

        var = StringVar(value = format_value(init_value, spec.units))
        self.__setattr__(f"{property}_str", var)

        name_label = Label(frame, text = spec.label + ":", anchor = "w", font = self._slider_font)
        name_label.grid(row = row, column = 0, sticky = "w", padx = (0, 6))

        value_label = Label(frame, textvariable = var, anchor = "e", width = 13, font = self._slider_font)
        value_label.grid(row = row, column = 1, sticky = "ew", padx = (0, 6))