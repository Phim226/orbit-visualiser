from tkinter import Tk, Frame, Scale, Label, StringVar, LabelFrame
from functools import partial
import numpy as np
from orbit_visualiser.ui import OrbitFigure
from orbit_visualiser.core import Orbit

class OrbitConfigurer():

    orbital_parameters: dict[str, str] = {
        "a" : "Semi-major axis",
        "b" : "Semi-minor axis",
        "ra": "Radius of apoapsis",
        "p" : "Semi-parameter"
    }

    def __init__(self, root: Tk, config_frame_placement: tuple[str], orbit_fig: OrbitFigure, orbit: Orbit):
        self._root = root

        self._display_frame = LabelFrame(root, text = "Orbital parameters")


        self._orbit_fig = orbit_fig
        self._orbit = orbit

        self._config_frame = Frame(root)
        self._config_frame.pack(side = config_frame_placement[0], anchor = config_frame_placement[1], padx = 8, pady = 6)

    def build(self) -> None:
        self._build_eccentricity_slider()
        self._build_periapsis_slider()

        self._display_frame.pack(side = "top", anchor = "nw", pady = (10, 0))

        for i, parameters in enumerate(list(self.orbital_parameters.items())):
            self._build_display(parameters[0], parameters[1], i)

    def _build_eccentricity_slider(self) -> None:
        self._e_slider = Scale(self._root, from_ = 0, to = 2, resolution = 0.01, length  = 150, orient = "horizontal",
                              command = partial(self._update_value, "e"), label = "Eccentricity", font=("Segoe UI", 9))
        self._e_slider.set(self._orbit.e)
        self._e_slider.pack(side = "top", anchor = "nw")

    def _build_periapsis_slider(self) -> None:
        self._rp_slider = Scale(self._root, from_ = 0, to = 10, resolution = 0.01, length  = 150, orient = "horizontal",
                              command = partial(self._update_value, "rp"), label = "Radius of periapsis", font=("Segoe UI", 9))
        self._rp_slider.set(self._orbit.rp)
        self._rp_slider.pack(side = "top", anchor = "nw")

    def _update_value(self, parameter: str, new_val: str) -> None:
        setattr(self._orbit, parameter, new_val)
        self._orbit_fig.redraw_orbit()

        for param, display_str in list(self.orbital_parameters.items()):
            self._update_display(param, display_str)

    def _build_display(self, parameter: str, display_str: str, row: int) -> None:
        var = StringVar(value = self._format_display_value(getattr(self._orbit, parameter)))
        self.__setattr__(f"_{parameter}_str", var)

        name_label = Label(self._display_frame, text = display_str + ":", anchor = "w", font=("Segoe UI", 9))
        name_label.grid(row = row, column = 0, sticky = "w", padx = (0, 6))

        value_label = Label(self._display_frame, textvariable = var, anchor = "e", width = 5, font=("Segoe UI", 9))
        value_label.grid(row = row, column = 1, sticky = "e", padx = (0, 6))

    def _update_display(self, parameter: str, display_str: str) -> None:
        self.__getattribute__(
            f"_{parameter}_str"
        ).set(self._format_display_value(getattr(self._orbit, parameter)))

    def _format_display_value(self, value: float) -> str:
        if value == np.inf:
            return "∞"

        return f"{value:6.2f}"