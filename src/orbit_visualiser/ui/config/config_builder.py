from tkinter import Tk, Frame, Scale, Label, StringVar, LabelFrame, Button, Entry, DoubleVar
from tkinter.ttk import Separator
from typing import Callable, TypeVar, Generic
from functools import partial
from dataclasses import dataclass
import numpy as np
from orbit_visualiser.ui import OrbitFigure
from orbit_visualiser.core import Orbit, Satellite, CentralBody

# TODO: Give option to show parameters on the plot (arrows/lines for vectors and distances etc).
# TODO: Split into variables, options and properties builders.
# TODO: Manage geometry of display options using rows/columns.


# TODO: VariableSpec can be subclass of PropertySpec
@dataclass(frozen = True)
class VariableSpec():
    label: str
    obj: Orbit | Satellite | CentralBody
    units: str | None
    init_value: float
    slider_lims: tuple[int]
    decimal_places: int
    entry_pos: tuple[int]

T = TypeVar("T")

@dataclass(frozen = True)
class PropertySpec(Generic[T]):
    label: str
    obj: Orbit | Satellite
    units: str | None
    getter: Callable[[T], float | bool]

class OrbitConfigBuilder():

    _title_font = ("Orbitron", 16, "bold")
    _subtitle_font = ("Orbitron", 11, "normal")
    _slider_font = ("Fira Mono", 9, "normal")

    orbital_properties: dict[str, PropertySpec[Orbit]] = {
        "orbit_type" : PropertySpec("Orbit type", None, lambda orbit: orbit.orbit_type),
        "semi_major_axis" : PropertySpec("Semi-major axis", "km", lambda orbit: orbit.a),
        "semi_minor_axis" : PropertySpec("Semi-minor axis", "km", lambda orbit: orbit.b),
        "radius_apoapsis": PropertySpec("Radius of apoapsis", "km", lambda orbit: orbit.ra),
        "semi_parameter" : PropertySpec("Semi-parameter", "km", lambda orbit: orbit.p),
        "asymptote_anomal" : PropertySpec("Anomaly of asymptote", "°", lambda orbit: np.degrees(orbit.t_asymp)),
        "turn_angle" : PropertySpec("Turning angle", "°", lambda orbit: np.degrees(orbit.turn_angle)),
        "aim_radius" : PropertySpec("Aiming radius", "km", lambda orbit: orbit.aim_rad)
    }

    satellite_properties: dict[str, PropertySpec[Satellite]] = {
        "radius" : PropertySpec("Radius", "km", lambda sat: sat.r),
        "x_pos" : PropertySpec("Perifocal x position", "km", lambda sat: sat.pos_pf[0]),
        "y_pos" : PropertySpec("Perifocal y position", "km", lambda sat: sat.pos_pf[1]),
        "period" : PropertySpec("Orbital period", "s", lambda sat: sat.period),
        "mean_motion" : PropertySpec("Mean motion", "°/s", lambda sat: np.degrees(sat.n)),
        "e_anomaly" : PropertySpec("Eccentric anomaly", "°", lambda sat: np.degrees(sat.e_anomaly)),
        "m_anomaly" : PropertySpec("Mean anomaly", "°", lambda sat: np.degrees(sat.m_anomaly)),
        "time_periapsis" : PropertySpec("Time since periapsis", "s", lambda sat: sat.t_p),
        "ang_momentum" : PropertySpec("Angular momentum", "km²/s", lambda sat: sat.h),
        "velocity" : PropertySpec("Velocity", "km/s", lambda sat: sat.v),
        "azim_velocity" : PropertySpec("Azimuthal velocity", "km/s", lambda sat: sat.v_azim),
        "radial_velocity" : PropertySpec("Radial velocity", "km/s", lambda sat: sat.v_radial),
        "esc_velocity" : PropertySpec("Escape velocity", "km/s", lambda sat: sat.v_esc),
        "excess_velocity" : PropertySpec("Excess velocity", "km/s", lambda sat: sat.v_inf),
        "flight_angle" : PropertySpec("Flight angle", "°", lambda sat: np.degrees(sat.gam)),
        "spec_energy" : PropertySpec("Specific energy", "km²/s²", lambda sat: sat.eps),
        "char_energy" : PropertySpec("Characteristic energy", "km²/s²", lambda sat: sat.c3)
    }

    properties: dict[str, PropertySpec] = orbital_properties | satellite_properties

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

        self._e_specs: VariableSpec = VariableSpec(
            "Eccentricity",
            orbit,
            None,
            orbit.e,
            (0, 5),
            3,
            (85, 4)
        )
        self._rp_specs: VariableSpec = VariableSpec(
            "Radius of periapsis",
            orbit,
            "km",
            orbit.rp,
            (central_body.r + 1, 200_000),
            0,
            (160, 4)
        )
        self._mu_specs: VariableSpec = VariableSpec(
            "Gravitational parameter",
            central_body,
            "km³/s²",
            central_body.mu,
            (1, 1_000_000),
            0,
            (198, 4)
        )
        self._nu_specs: VariableSpec = VariableSpec(
            "True anomaly",
            satellite,
            "°",
            satellite.nu,
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

        self._property_specs_by_object: dict[Orbit | Satellite, dict] = {
            orbit: self.orbital_properties,
            satellite : self.satellite_properties
        }

        self._config_frame = Frame(root)
        self._config_frame.pack(
            side = config_frame_placement[0],
            anchor = config_frame_placement[1],
            padx = 8, pady = 6
        )

    @property
    def variable_specs(self) -> dict[str, VariableSpec]:
        return self._variable_specs

    @property
    def property_specs_by_object(self) -> dict[Orbit | Satellite, dict]:
        return self._property_specs_by_object

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
            parameter: str,
            param_props: VariableSpec,
            validate_input: Callable,
            update_value: Callable
    ) -> tuple[Scale, Entry]:
        obj = param_props.obj
        units = param_props.units

        frame = Frame(root, width = 265, height = 60)

        slider = self._build_slider(
            frame,
            parameter,
            obj,
            f"{param_props.name}{"" if units is None else f" ({units})"} = ",
            param_props.slider_lims,
            1/10**param_props.decimal_places,
            update_value
        )

        entry = Entry(frame, width = 10)
        entry.insert(0, f"{param_props.init_value: 0.{param_props.decimal_places}f}".strip())
        entry.bind("<Return>", partial(validate_input, parameter, obj))
        x, y = param_props.entry_pos
        entry.place(x = x, y = y)

        frame.pack(side = "top", anchor = "nw", pady = 2)

        return slider, entry

    def _build_slider(
            self,
            root: Frame,
            parameter: str,
            source_object: Orbit | Satellite,
            label: str,
            lims: tuple[int],
            res: float,
            update_value: Callable
    ) -> Scale:
        slider_var: DoubleVar = DoubleVar()
        self.__setattr__(f"{parameter}_var", slider_var)

        slider_name = f"_{parameter}_slider"
        self.__setattr__(
            slider_name,
            Scale(root, from_ = lims[0], to = lims[1], resolution = res, length = 260,
                  orient = "horizontal", variable = slider_var,
                  command = partial(update_value, parameter, source_object, "slider"),
                  label = label, font = self._slider_font)
        )

        init_value: float = (round(np.degrees(getattr(source_object, parameter)), 2) if parameter == "nu"
                             else getattr(source_object, parameter))
        slider_var.set(init_value)

        slider: Scale = self.__getattribute__(slider_name)
        slider.place(x = 0, y = 0, anchor = "nw")
        return slider

    def _build_properties_frame(self, format_value: Callable) -> None:
        props_frame = Frame(self._config_frame, padx = 2)
        self._properties_frame = props_frame

        self._build_separator(props_frame, "Properties")
        orbital_props_frame = LabelFrame(
            props_frame, bd = 1, relief = "sunken", text = "Orbit", font = self._subtitle_font
        )
        self._populate_properties(orbital_props_frame, self.orbital_properties, self._orbit, format_value)
        orbital_props_frame.pack(side = "top", anchor = "nw", pady = (2, 0))

        sat_props_frame = LabelFrame(
            props_frame, bd = 1, relief = "sunken", text = "Satellite",
            font = self._subtitle_font, width = 244
        )
        self._populate_properties(sat_props_frame, self.satellite_properties, self._sat, format_value)
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

    def _build_separator(self, root: Frame, text: str) -> None:
        frame = Frame(root)
        frame.pack(side = "top", fill = "x", pady = 4)

        Label(frame, text = text, font = self._title_font).pack(side = "left", padx = (0, 6))
        Frame(frame, height = 2, bd = 1, relief = "sunken").pack(side = "left", fill = "x", expand = True)