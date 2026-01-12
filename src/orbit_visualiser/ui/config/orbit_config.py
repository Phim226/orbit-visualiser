from tkinter import Tk, Frame, Scale, Label, StringVar, LabelFrame
from functools import partial
import numpy as np
from orbit_visualiser.ui import OrbitFigure
from orbit_visualiser.core import Orbit, Satellite, CentralBody

class OrbitConfigurer():

    orbital_parameters: dict[str, tuple[str]] = {
        "a" : ("Semi-major axis", "km"),
        "b" : ("Semi-minor axis", "km"),
        "ra": ("Radius of apoapsis", "km"),
        "p" : ("Semi-parameter", "km"),
        "t_asymp" : ("Asymptote anomaly", "°"),
        "turn_angle" : ("Turning angle", "°"),
        "aim_rad" : ("Aiming radius", "km")
    }

    satellite_parameters: dict[str, tuple[str]] = {

    }

    parameters: dict[str, tuple[str]] = orbital_parameters | satellite_parameters

    def __init__(self, root: Tk, config_frame_placement: tuple[str], orbit_fig: OrbitFigure, orbit: Orbit, central_body: CentralBody, satellite: Satellite):
        self._root = root

        self._orbit_fig = orbit_fig
        self._orbit = orbit
        self._central_body = central_body
        self._sat = satellite

        self._config_frame = Frame(root)
        self._config_frame.pack(side = config_frame_placement[0], anchor = config_frame_placement[1], padx = 8, pady = 6)

        self._slider_frame = Frame(self._config_frame)
        self._display_frame = Frame(self._config_frame)

    def build(self) -> None:
        self._e_slider = self._build_slider("e", self._orbit, "Eccentricity", 2, 0.01)
        self._rp_slider = self._build_slider("rp", self._orbit, "Radius of periapsis (km)", 10_000)
        self._mu_slider = self._build_slider("mu", self._central_body, "Gravitational parameter (km³/s²)", 1_000_000)
        self._nu_slider = self._build_slider("nu", self._sat, "True anomaly (°)", 360)

        self._slider_frame.pack(side = "top", anchor = "nw", pady = (10, 0))
        self._display_frame.pack(side = "top", anchor = "nw", pady = (10, 0))

        for i, parameters in enumerate(list(self.parameters.items())):
            parameter_info = parameters[1]
            self._build_display(parameters[0], parameter_info[0], parameter_info[1], i)

    def _build_slider(self, parameter: str, source_object: Orbit | Satellite, label: str, upper_lim: int, res: float = 1) -> Scale:
        slider_name = f"_{parameter}_slider"
        self.__setattr__(
            slider_name,
            Scale(self._slider_frame, to = upper_lim, resolution = res, length = 175, orient = "horizontal",
                  command = partial(self._update_value, parameter, source_object), label = label, font = ("Segoe UI", 9))
        )
        slider: Scale = self.__getattribute__(slider_name)
        slider.set(getattr(source_object, parameter))
        slider.pack(side = "top", anchor = "nw")
        return slider

    def _update_value(self, parameter: str, source_object: Orbit | Satellite, new_val: str) -> None:
        new_val = float(new_val)
        setattr(source_object, parameter, new_val)
        self._orbit_fig.redraw_orbit()

        if parameter == "e":
            if new_val > 1:
                t_asymp = self._orbit.t_asymp
                self._nu_slider.configure(from_ = -t_asymp, to = t_asymp)
            else:
                self._nu_slider.configure(from_ = 0, to = 360)

        for param in self.parameters:
            self._update_display(param)

    def _build_display(self, parameter: str, display_str: str, units: str, row: int) -> None:
        var = StringVar(value = self._format_display_value(getattr(self._orbit, parameter), units))
        self.__setattr__(f"_{parameter}_str", var)

        name_label = Label(self._display_frame, text = display_str + ":", anchor = "w", font=("Segoe UI", 9))
        name_label.grid(row = row, column = 0, sticky = "w", padx = (0, 6))

        value_label = Label(self._display_frame, textvariable = var, anchor = "e", width = 10, font=("Segoe UI", 9))
        value_label.grid(row = row, column = 1, sticky = "e", padx = (0, 6))

    def _update_display(self, parameter: str) -> None:
        self.__getattribute__(
            f"_{parameter}_str"
        ).set(self._format_display_value(getattr(self._orbit, parameter), self.parameters[parameter][1]))

    def _format_display_value(self, value: float, units: str) -> str:
        if np.isinf(value):
            return f"∞ {units}"

        elif np.isnan(value):
            return "n/a"

        elif units == "km":
            return f"{value:6.0f} km"

        elif units == "°":
            return f"{value:6.2f} °"