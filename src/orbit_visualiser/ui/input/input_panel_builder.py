from tkinter import Frame, Scale, LabelFrame, Button, Entry, DoubleVar, Label
from typing import Callable
from functools import partial
import numpy as np
from orbit_visualiser.ui.common.builder import Builder
from orbit_visualiser.ui.common.specs import VariableSpec
from orbit_visualiser.ui.common.presets import initial_config
from orbit_visualiser.ui.data_access import OrbitDataAccess


class InputBuilder(Builder):

    def __init__(
            self,
            input_frame: Frame,
            oda: OrbitDataAccess
    ):
        self._input_frame = input_frame

        self._oda = oda

        self._e_specs: VariableSpec = VariableSpec(
            "Eccentricity",
            None,
            lambda sat: sat.orbit.eccentricity,
            initial_config.eccentricity,
            (0, 5),
            3,
            "normal"
        )
        self._rp_specs: VariableSpec = VariableSpec(
            "Radius of periapsis",
            "km",
            lambda sat: sat.orbit.radius_of_periapsis,
            initial_config.radius_of_periapsis,
            (initial_config.radius + 1, 200_000),
            0,
            "normal"
        )
        self._mu_specs: VariableSpec = VariableSpec(
            "Gravitational parameter",
            "km³/s²",
            lambda sat: sat.central_body.mu,
            initial_config.gravitational_parameter,
            (1, 1_000_000),
            0,
            "normal"
        )
        self._nu_specs: VariableSpec = VariableSpec(
            "True anomaly",
            "°",
            lambda sat: np.degrees(sat.true_anomaly),
            np.degrees(initial_config.true_anomaly),
            (0, 360),
            2,
            "normal"
        )
        self._raan_specs: VariableSpec = VariableSpec(
            "Right ascension of the ascending node",
            "°",
            lambda sat: np.degrees(sat.orbit.right_ascen_of_ascend_node),
            np.degrees(initial_config.right_ascension_of_the_ascending_node),
            (0, 360),
            2,
            "disabled"
        )
        self._i_specs: VariableSpec = VariableSpec(
            "Inclination",
            "°",
            lambda sat: np.degrees(sat.orbit.inclination),
            np.degrees(initial_config.inclination),
            (0, 180),
            2,
            "normal"
        )
        self._omega_specs: VariableSpec = VariableSpec(
            "Argument of periapsis",
            "°",
            lambda sat: np.degrees(sat.orbit.argument_of_periapsis),
            np.degrees(initial_config.argument_of_periapsis),
            (0, 360),
            2,
            "disabled"
        )

        self._variable_specs: dict[str, VariableSpec] = {
            "e" : self._e_specs,
            "rp" : self._rp_specs,
            "nu" : self._nu_specs,
            "raan" : self._raan_specs,
            "i" : self._i_specs,
            "omega" : self._omega_specs,
            "mu" : self._mu_specs
        }

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

    @property
    def raan_slider(self) -> Scale:
        return self._raan_slider

    @property
    def raan_entry(self) -> Entry:
        return self._raan_entry

    @property
    def i_slider(self) -> Scale:
        return self._i_slider

    @property
    def i_entry(self) -> Entry:
        return self._i_entry

    @property
    def omega_slider(self) -> Scale:
        return self._omega_slider

    @property
    def omega_entry(self) -> Entry:
        return self._omega_entry

    def build_input_frame(
            self,
            reset: Callable,
            validate_input: Callable,
            slider_changed: Callable
    ) -> None:
        var_frame = Frame(self._input_frame)
        self._variables_frame = var_frame

        self._build_separator(var_frame, "Variables")

        # Build orbital geometry frame
        orbital_geom_frame = LabelFrame(
            var_frame, bd = 2, relief = "sunken", text = "Orbital geometry", font = self._subtitle_font
        )
        self._e_slider, self._e_entry = self._build_input_frame(
            orbital_geom_frame, "e", self._variable_specs["e"], validate_input, slider_changed
        )
        self._rp_slider, self._rp_entry = self._build_input_frame(
            orbital_geom_frame, "rp", self._variable_specs["rp"], validate_input, slider_changed
        )
        self._raan_slider, self._raan_entry = self._build_input_frame(
            orbital_geom_frame, "raan", self._variable_specs["raan"], validate_input, slider_changed
        )
        self._i_slider, self._i_entry = self._build_input_frame(
            orbital_geom_frame, "i", self._variable_specs["i"], validate_input, slider_changed
        )
        self._omega_slider, self._omega_entry = self._build_input_frame(
            orbital_geom_frame, "omega", self._variable_specs["omega"], validate_input, slider_changed
        )
        orbital_geom_frame.pack(side = "top", anchor = "nw", pady = (4, 0))

        # Building central body frame
        central_body_frame = LabelFrame(
            var_frame, bd = 2, relief = "sunken", text = "Central body", font = self._subtitle_font
        )
        self._mu_slider, self._mu_entry = self._build_input_frame(
            central_body_frame, "mu", self._variable_specs["mu"], validate_input, slider_changed
        )
        central_body_frame.pack(side = "top", anchor = "nw", pady = (4, 0))

        # Build satellite frame
        sat_frame = LabelFrame(
            var_frame, bd = 2, relief = "sunken", text = "Satellite", font = self._subtitle_font
        )
        self._nu_slider, self._nu_entry = self._build_input_frame(
            sat_frame, "nu", self._variable_specs["nu"], validate_input, slider_changed
        )
        sat_frame.pack(side = "top", anchor = "nw", pady = (4, 0))

        # Build reset button
        reset_button = Button(var_frame, text = "Reset", command = reset)
        reset_button.pack(side = "top", anchor = "nw", pady = (4, 0))

        var_frame.pack(side = "top", anchor = "nw", pady = (2, 0))

    def _build_input_frame(
            self,
            root: Frame,
            variable: str,
            spec: VariableSpec,
            validate_input: Callable,
            slider_changed: Callable
    ) -> tuple[Scale, Entry]:
        frame = Frame(root, width = 290, height = 75, relief = "groove", bd = 1)

        units = spec.units
        label = Label(frame, text = f"{spec.label}{"" if units is None else f" ({units})"}:")
        label.place(x = 5, y = 0)

        slider = self._build_slider(
            frame,
            variable,
            spec,
            slider_changed
        )

        entry = Entry(frame, width = 10)
        entry.insert(0, f"{spec.getter(self._oda.satellite): 0.{spec.decimal_places}f}".strip())
        entry.configure(state = spec.init_state)
        entry.bind("<Return>", partial(validate_input, variable))
        entry.place(x = 5, y = 20)

        frame.pack(side = "top", anchor = "nw", pady = 2)

        return slider, entry

    def _build_slider(
            self,
            root: Frame,
            variable: str,
            spec: VariableSpec,
            slider_changed: Callable
    ) -> Scale:
        slider_var: DoubleVar = DoubleVar()
        self.__setattr__(f"{variable}_var", slider_var)

        slider_name = f"_{variable}_slider"
        lims = spec.slider_lims
        self.__setattr__(
            slider_name,
            Scale(root, from_ = lims[0], to = lims[1], resolution = 1/10**spec.decimal_places, length = 275,
                  orient = "horizontal", variable = slider_var,
                  command = partial(slider_changed, variable, "slider"),
                  tickinterval = 0, showvalue = 0,
                  state = spec.init_state
                  )
        )

        slider_var.set(spec.getter(self._oda.satellite))

        slider: Scale = self.__getattribute__(slider_name)
        slider.place(x = 5, y = 45, anchor = "nw")
        return slider