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
        self._build_slider("e", "Eccentricity", 2, 0.01)
        self._build_slider("rp", "Radius of periapsis (km)", 10000)

        self._display_frame.pack(side = "top", anchor = "nw", pady = (10, 0))

        for i, parameters in enumerate(list(self.orbital_parameters.items())):
            self._build_display(parameters[0], parameters[1], i)

    def _build_slider(self, parameter: str, label: str, upper_lim: int, res: float = 1) -> None:
        slider_name = f"_{parameter}_slider"
        self.__setattr__(
            slider_name,
            Scale(self._root, to = upper_lim, resolution = res, length = 150, orient = "horizontal",
                  command = partial(self._update_value, parameter), label = label, font = ("Segoe UI", 9))
        )
        slider: Scale = self.__getattribute__(slider_name)
        slider.set(getattr(self._orbit, parameter))
        slider.pack(side = "top", anchor = "nw")

    def _update_value(self, parameter: str, new_val: str) -> None:
        setattr(self._orbit, parameter, float(new_val))
        self._orbit_fig.redraw_orbit()

        for param in self.orbital_parameters:
            self._update_display(param)

    def _build_display(self, parameter: str, display_str: str, row: int) -> None:
        var = StringVar(value = f"{self._format_display_value(getattr(self._orbit, parameter))} km")
        self.__setattr__(f"_{parameter}_str", var)

        name_label = Label(self._display_frame, text = display_str + ":", anchor = "w", font=("Segoe UI", 9))
        name_label.grid(row = row, column = 0, sticky = "w", padx = (0, 6))

        value_label = Label(self._display_frame, textvariable = var, anchor = "e", width = 8, font=("Segoe UI", 9))
        value_label.grid(row = row, column = 1, sticky = "e", padx = (0, 6))

    def _update_display(self, parameter: str) -> None:
        self.__getattribute__(f"_{parameter}_str").set(f"{self._format_display_value(getattr(self._orbit, parameter))} km")

    def _format_display_value(self, value: float) -> str:
        if value == np.inf:
            return "∞"

        return f"{value:6.0f}"