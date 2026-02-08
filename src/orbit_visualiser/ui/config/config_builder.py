from tkinter import Tk, Frame, Scale, Label, StringVar, LabelFrame, Button, Entry, DoubleVar, Checkbutton, IntVar
from tkinter.ttk import Separator
from typing import Callable
from functools import partial
import numpy as np
from orbit_visualiser.core import Orbit, Satellite, CentralBody
from orbit_visualiser.ui.config.display_panel.display_panel_builder import DisplayBuilder
from orbit_visualiser.ui.config.properties_panel.properties_panel_builder import PropertiesBuilder
from orbit_visualiser.ui.config.variables_panel.variables_panel_builder import VariablesBuilder
from orbit_visualiser.ui.common.builder import Builder
from orbit_visualiser.ui.common.specs import PropertySpec, VariableSpec

# TODO: Give option to show parameters on the plot (arrows/lines for vectors and distances etc).
# TODO: Split into variables, options and properties builders.
# TODO: Manage geometry of display options using rows/columns.

class OrbitConfigBuilderTest(Builder):

    def __init__(
            self, root: Tk,
            config_frame_placement: tuple[str],
            orbit: Orbit,
            central_body: CentralBody,
            satellite: Satellite
    ):
        self._root = root

        self._orbit = orbit
        self._central_body = central_body
        self._sat = satellite

        self._config_frame = Frame(root)
        self._config_frame.pack(
            side = config_frame_placement[0],
            anchor = config_frame_placement[1],
            padx = 8, pady = 6
        )

        self._options_frame = Frame(self._config_frame, padx = 2)
        self._options_frame.pack(side = "left", anchor = "n", pady = (2, 0))

        self._variables_builder = VariablesBuilder(self._options_frame, orbit, central_body, satellite)
        self._display_builder = DisplayBuilder(self._options_frame)
        self._properties_builder = PropertiesBuilder(self._config_frame, orbit, satellite)

    @property
    def variables_builder(self) -> VariablesBuilder:
        return self._variables_builder

    @property
    def display_builder(self) -> DisplayBuilder:
        return self._display_builder

    @property
    def properties_builder(self) -> PropertiesBuilder:
        return self._properties_builder

    def build(
            self,
            reset: Callable,
            validate_input: Callable,
            slider_changed: Callable,
            format_value: Callable
    ) -> None:


        self._variables_builder.build_variables_frame(reset, validate_input, update_value)
        self._display_builder.build_display_options_frame()


        sep = Separator(self._config_frame, orient = "vertical")
        sep.pack(side = "left", fill = "y", padx = 6, expand = True)

        self._properties_builder.build_properties_frame(format_value)


class OrbitConfigBuilder(Builder):

    def __init__(
            self, root: Tk,
            config_frame_placement: tuple[str],
            orbit: Orbit,
            central_body: CentralBody,
            satellite: Satellite
    ):
        self._root = root

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

        self._e_specs: VariableSpec = VariableSpec(
            "Eccentricity",
            orbit,
            None,
            lambda orbit: orbit.e,
            orbit.e,
            (0, 5),
            3,
            (85, 4)
        )
        self._rp_specs: VariableSpec = VariableSpec(
            "Radius of periapsis",
            orbit,
            "km",
            lambda orbit: orbit.rp,
            orbit.rp,
            (central_body.r + 1, 200_000),
            0,
            (160, 4)
        )
        self._mu_specs: VariableSpec = VariableSpec(
            "Gravitational parameter",
            central_body,
            "km³/s²",
            lambda central_body: central_body.mu,
            central_body.mu,
            (1, 1_000_000),
            0,
            (198, 4)
        )
        self._nu_specs: VariableSpec = VariableSpec(
            "True anomaly",
            satellite,
            "°",
            lambda sat: np.degrees(sat.nu),
            np.degrees(satellite.nu),
            (0, 360),
            2,
            (115, 4)
        )

        self._variable_specs: dict[str, VariableSpec] = {
            "e" : self._e_specs,
            "rp" : self._rp_specs,
            "mu" : self._mu_specs,
            "nu" : self._nu_specs
        }

        self._config_frame = Frame(root)
        self._config_frame.pack(
            side = config_frame_placement[0],
            anchor = config_frame_placement[1],
            padx = 8, pady = 6
        )

    @property
    def property_specs(self) -> dict[str, PropertySpec]:
        return self._property_specs

    @property
    def variable_specs(self) -> dict[str, VariableSpec]:
        return self._variable_specs

    @property
    def e_slider(self) -> Scale:
        return self._e_slider

    @property
    def e_entry(self) -> Entry:
        return self._e_entry

    @property
    def rp_slider(self) -> Scale:
        return self._rp_slider

    @property
    def rp_entry(self) -> Entry:
        return self._rp_entry

    @property
    def mu_slider(self) -> Scale:
        return self._mu_slider

    @property
    def mu_entry(self) -> Entry:
        return self._mu_entry

    @property
    def nu_slider(self) -> Scale:
        return self._nu_slider

    @property
    def nu_entry(self) -> Entry:
        return self._nu_entry

    def build(
            self,
            reset: Callable,
            validate_input: Callable,
            update_value: Callable,
            format_value: Callable
    ) -> None:
        self._options_frame = Frame(self._config_frame, padx = 2)
        self._options_frame.pack(side = "left", anchor = "n", pady = (2, 0))

        self._build_variables_frame(reset, validate_input, update_value)
        self._build_display_options_frame()


        sep = Separator(self._config_frame, orient = "vertical")
        sep.pack(side = "left", fill = "y", padx = 6, expand = True)

        self._build_properties_frame(format_value)

    def _build_variables_frame(
            self,
            reset: Callable,
            validate_input: Callable,
            update_value: Callable
    ) -> None:
        var_frame = Frame(self._options_frame)
        self._variables_frame = var_frame

        self._build_separator(var_frame, "Variables")

        # Build orbital geometry frame
        orbital_geom_frame = LabelFrame(
            var_frame, bd = 1, relief = "sunken", text = "Orbital geometry", font = self._subtitle_font
        )
        self._e_slider, self._e_entry = self._build_input_frame(
            orbital_geom_frame, "e", self._variable_specs["e"], validate_input, update_value
        )
        self._rp_slider, self._rp_entry = self._build_input_frame(
            orbital_geom_frame, "rp", self._variable_specs["rp"], validate_input, update_value
        )
        orbital_geom_frame.pack(side = "top", anchor = "nw", pady = (4, 0))

        # Building central body frame
        central_body_frame = LabelFrame(
            var_frame, bd = 1, relief = "sunken", text = "Central body", font = self._subtitle_font
        )
        self._mu_slider, self._mu_entry = self._build_input_frame(
            central_body_frame, "mu", self._variable_specs["mu"], validate_input, update_value
        )
        central_body_frame.pack(side = "top", anchor = "nw", pady = (4, 0))

        # Build satellite frame
        sat_frame = LabelFrame(
            var_frame, bd = 1, relief = "sunken", text = "Satellite", font = self._subtitle_font
        )
        self._nu_slider, self._nu_entry = self._build_input_frame(
            sat_frame, "nu", self._variable_specs["nu"], validate_input, update_value
        )
        sat_frame.pack(side = "top", anchor = "nw", pady = (4, 0))

        # Build reset button
        reset_button = Button(var_frame, text = "Reset", command = reset)
        reset_button.pack(side = "top", anchor = "nw", pady = (4, 0))

        var_frame.pack(side = "top", anchor = "ne", pady = (2, 0))

    def _build_input_frame(
            self,
            root: Frame,
            variable: str,
            spec: VariableSpec,
            validate_input: Callable,
            update_value: Callable
    ) -> tuple[Scale, Entry]:
        obj = spec.obj
        units = spec.units

        frame = Frame(root, width = 265, height = 60)

        slider = self._build_slider(
            frame,
            variable,
            spec,
            update_value
        )

        entry = Entry(frame, width = 10)
        entry.insert(0, f"{spec.getter(obj): 0.{spec.decimal_places}f}".strip())
        entry.bind("<Return>", partial(validate_input, variable, obj))
        x, y = spec.entry_pos
        entry.place(x = x, y = y)

        frame.pack(side = "top", anchor = "n", pady = 2)

        return slider, entry

    def _build_slider(
            self,
            root: Frame,
            variable: str,
            spec: VariableSpec,
            update_value: Callable
    ) -> Scale:
        slider_var: DoubleVar = DoubleVar()
        self.__setattr__(f"{variable}_var", slider_var)

        slider_name = f"_{variable}_slider"
        lims = spec.slider_lims
        obj = spec.obj
        units = spec.units
        label = f"{spec.label}{"" if units is None else f" ({units})"} = "
        self.__setattr__(
            slider_name,
            Scale(root, from_ = lims[0], to = lims[1], resolution = 1/10**spec.decimal_places, length = 260,
                  orient = "horizontal", variable = slider_var,
                  command = partial(update_value, variable, obj, "slider"),
                  label = label, font = self._slider_font)
        )

        slider_var.set(spec.getter(obj))

        slider: Scale = self.__getattribute__(slider_name)
        slider.place(x = 0, y = 0, anchor = "nw")
        return slider

    def _build_display_options_frame(self) -> None:
        display_options_frame: Frame = Frame(self._options_frame)
        self._display_options_frame = display_options_frame

        self._build_separator(display_options_frame, "Display options")

        # Build orbit options frame
        orbit_options_frame = LabelFrame(
            display_options_frame, bd = 1, relief = "sunken", text = "Orbit", font = self._subtitle_font
        )
        self._rp_display_var: IntVar = IntVar()
        rp_display_check: Checkbutton = Checkbutton(orbit_options_frame, text = "Periapsis", variable = self._rp_display_var)
        rp_display_check.pack(side = "top", anchor = "nw")
        orbit_options_frame.pack(side = "top", anchor = "nw", pady = (4, 0))

        # Build central body options frame
        central_body_options_frame = LabelFrame(
            display_options_frame, bd = 1, relief = "sunken", text = "Central body", font = self._subtitle_font
        )
        self._central_body_display_var: IntVar = IntVar(value = 1)
        central_body_display_check: Checkbutton = Checkbutton(central_body_options_frame, text = "Central body", variable = self._central_body_display_var)
        central_body_display_check.pack(side = "top", anchor = "nw")
        central_body_options_frame.pack(side = "top", anchor = "nw", pady = (4, 0))

        # Build satellite options frame
        satellite_options_frame = LabelFrame(
            display_options_frame, bd = 1, relief = "sunken", text = "Satellite", font = self._subtitle_font
        )
        self._radius_display_var: IntVar = IntVar()
        radius_display_check: Checkbutton = Checkbutton(satellite_options_frame, text = "Radius", variable = self._radius_display_var)
        radius_display_check.pack(side = "top", anchor = "nw")
        satellite_options_frame.pack(side = "top", anchor = "nw", pady = (4, 0))

        display_options_frame.pack(side = "top", anchor = "nw", pady = (2, 0), fill = "x", expand = True)

    def _build_properties_frame(self, format_value: Callable) -> None:
        props_frame = Frame(self._config_frame, padx = 2)
        self._properties_frame = props_frame

        self._build_separator(props_frame, "Properties")
        orbital_props_frame = LabelFrame(
            props_frame, bd = 1, relief = "sunken", text = "Orbit", font = self._subtitle_font
        )
        self._populate_properties(orbital_props_frame, self._orbital_properties, self._orbit, format_value)
        orbital_props_frame.pack(side = "top", anchor = "nw", pady = (2, 0))

        sat_props_frame = LabelFrame(
            props_frame, bd = 1, relief = "sunken", text = "Satellite",
            font = self._subtitle_font, width = 244
        )
        self._populate_properties(sat_props_frame, self._satellite_properties, self._sat, format_value)
        sat_props_frame.pack(side = "top", anchor = "nw", pady = (2, 0), fill = "x")

        props_frame.pack(side = "top", anchor = "n", pady = (2, 0))

    def _populate_properties(
            self,
            frame: LabelFrame,
            properties: dict[str, PropertySpec],
            source_object: Orbit | Satellite,
            format_value: Callable
    ) -> None:
        for i, (key, spec) in enumerate(properties.items()):
            self._build_property_row(
                frame, key, source_object, spec, i, format_value
        )

        # Setting the weight allows the grid manager to stretch labels in _build_display into available space.
        frame.grid_columnconfigure(0, weight = 0)
        frame.grid_columnconfigure(1, weight = 1)

    def _build_property_row(
            self,
            frame: LabelFrame,
            property: str,
            source_object: Orbit | Satellite,
            spec: PropertySpec,
            row: int,
            format_value: Callable
    ) -> None:
        init_value = spec.getter(source_object)

        var = StringVar(value = format_value(init_value, spec.units))
        self.__setattr__(f"{property}_str", var)

        name_label = Label(frame, text = spec.label + ":", anchor = "w", font = self._slider_font)
        name_label.grid(row = row, column = 0, sticky = "w", padx = (0, 6))

        value_label = Label(frame, textvariable = var, anchor = "e", width = 13, font = self._slider_font)
        value_label.grid(row = row, column = 1, sticky = "ew", padx = (0, 6))