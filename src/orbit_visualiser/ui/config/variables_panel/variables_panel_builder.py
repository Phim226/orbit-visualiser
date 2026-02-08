from tkinter import Tk, Frame, Scale, LabelFrame, Button, Entry, DoubleVar
from typing import Callable
from functools import partial
import numpy as np
from orbit_visualiser.core import Orbit, Satellite, CentralBody
from orbit_visualiser.ui.common.builder import Builder
from orbit_visualiser.ui.common.specs import VariableSpec


class VariablesBuilder(Builder):


    def __init__(
            self,
            options_frame: Frame,
            orbit: Orbit,
            central_body: CentralBody,
            satellite: Satellite
    ):
        self._options_frame = options_frame

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

    def build_variables_frame(
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