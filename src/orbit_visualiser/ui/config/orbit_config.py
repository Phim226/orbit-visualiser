from tkinter import Tk, Frame, Scale, Label, StringVar, LabelFrame, Button
from tkinter.ttk import Separator
from functools import partial
import numpy as np
from orbit_visualiser.ui import OrbitFigure
from orbit_visualiser.core import Orbit, Satellite, CentralBody

# TODO: Display satellite parameters.
# TODO: Give option to show parameters on the plot (arrows/lines for vectors and distances etc).
class OrbitConfigurer():

    title_font = ("Orbitron", 16, "bold")
    subtitle_font = ("Orbitron", 11, "normal")
    slider_font = ("Fira Mono", 9, "normal")

    orbital_parameters: dict[str, tuple[str]] = {
        "a" : ("Semi-major axis", "km"),
        "b" : ("Semi-minor axis", "km"),
        "ra": ("Radius of apoapsis", "km"),
        "p" : ("Semi-parameter", "km"),
        "t_asymp" : ("Anomaly of asymptote", "°"),
        "turn_angle" : ("Turning angle", "°"),
        "aim_rad" : ("Aiming radius", "km")
    }

    satellite_parameters: dict[str, tuple[str]] = {
        "r" : ("Radius", "km"),
        "h" : ("Angular momentum", "km²/s"),
        "eps" : ("Mechanical energy", "km²/s²"),
        "v" : ("Velocity", "km/s"),
        "v_azim" : ("Azimuthal velocity", "km/s"),
        "v_radial" : ("Radial velocity", "km/s"),
        "v_esc" : ("Escape velocity", "km/s"),
        "v_inf" : ("Excess velocity", "km/s")
    }

    parameters: dict[str, tuple[str]] = orbital_parameters | satellite_parameters

    def __init__(self, root: Tk, config_frame_placement: tuple[str], orbit_fig: OrbitFigure, orbit: Orbit, central_body: CentralBody, satellite: Satellite):
        self._root = root

        self._orbit_fig = orbit_fig
        self._orbit = orbit
        self._central_body = central_body
        self._sat = satellite

        self._initial_state: dict[str, float] = {
            "e" : orbit.e,
            "rp" : orbit.rp,
            "mu" : central_body.mu,
            "nu" : satellite.nu
        }

        self._variable_objects: dict[str: Orbit | Satellite | CentralBody] = {
            "e" : orbit,
            "rp": orbit,
            "mu" : central_body,
            "nu" : satellite
        }

        self._parameter_objects: dict[Orbit | Satellite, dict] = {orbit: self.orbital_parameters, satellite : self.satellite_parameters}

        self._config_frame = Frame(root)
        self._config_frame.pack(side = config_frame_placement[0], anchor = config_frame_placement[1], padx = 8, pady = 6)

    def build(self) -> None:
        self._build_variables_frame()

        sep = Separator(self._config_frame, orient = "vertical")
        sep.pack(side = "left", fill = "y", padx = 6, expand = True)

        self._build_properties_frame()

    def _build_variables_frame(self) -> None:
        var_frame = Frame(self._config_frame, padx = 2)
        self._variables_frame = var_frame

        self._build_separator(var_frame, "Variables")
        orbital_geom_frame = LabelFrame(var_frame, bd = 1, relief = "sunken", text = "Orbital geometry", font = self.subtitle_font)
        self._e_slider = self._build_slider(orbital_geom_frame, "e", self._orbit, "Eccentricity", 2, res = 0.01)
        self._rp_slider = self._build_slider(orbital_geom_frame, "rp", self._orbit, "Radius of periapsis (km)", 100_000, lower_lim = self._central_body.r + 1)
        orbital_geom_frame.pack(side = "top", anchor = "nw", pady = (4, 0))

        attracting_body_frame = LabelFrame(var_frame, bd = 1, relief = "sunken", text = "Central body", font = self.subtitle_font)
        self._mu_slider = self._build_slider(attracting_body_frame, "mu", self._central_body, "Gravitational parameter (km³/s²)", 1_000_000, lower_lim = 1)
        attracting_body_frame.pack(side = "top", anchor = "nw", pady = (4, 0))

        sat_frame = LabelFrame(var_frame, bd = 1, relief = "sunken", text = "Satellite", font = self.subtitle_font)
        self._nu_slider = self._build_slider(sat_frame, "nu", self._sat, "True anomaly (°)", 360, res = 0.01)
        sat_frame.pack(side = "top", anchor = "nw", pady = (4, 0))

        reset_button = Button(var_frame, text = "Reset", command = self._reset_state)
        reset_button.pack(side = "top", anchor = "nw", pady = (4, 0))

        var_frame.pack(side = "left", anchor = "n", pady = (2, 0))

    def _reset_state(self) -> None:
        for name, value in list(self._initial_state.items()):
            self.__getattribute__(f"_{name}_slider").set(value)
            setattr(self._variable_objects[name], name, value)

        self._orbit.update_orbital_properties()
        self._orbit.update_orbit_type()
        self._sat.update_satellite_properties()
        self._orbit_fig.redraw_orbit()

    def _build_properties_frame(self) -> None:
        props_frame = Frame(self._config_frame, padx = 2)
        self._properties_frame = props_frame

        self._build_separator(props_frame, "Properties")
        orbital_props_frame = LabelFrame(props_frame, bd = 1, relief = "sunken", text = "Orbit", font = self.subtitle_font)
        self._populate_properties(orbital_props_frame, self.orbital_parameters, self._orbit)
        orbital_props_frame.pack(side = "top", anchor = "nw", pady = (2, 0))

        sat_props_frame = LabelFrame(props_frame, bd = 1, relief = "sunken", text = "Satellite", font = self.subtitle_font, width = 244)
        self._populate_properties(sat_props_frame, self.satellite_parameters, self._sat)
        sat_props_frame.pack(side = "top", anchor = "nw", pady = (2, 0), fill = "x")

        props_frame.pack(side = "top", anchor = "n", pady = (2, 0))

    def _populate_properties(self, frame: LabelFrame, parameters: dict[str, tuple[str]], source_object: Orbit | Satellite) -> None:
        for i, parameters in enumerate(list(parameters.items())):
            parameter_info = parameters[1]
            self._build_display(frame, parameters[0], source_object, parameter_info[0], parameter_info[1], i)

        # Setting the weight allows the grid manager to stretch labels in _build_display into available space.
        frame.grid_columnconfigure(0, weight=0)
        frame.grid_columnconfigure(1, weight=1)

    def _build_separator(self, root: Frame, text: str) -> None:
        frame = Frame(root)
        frame.pack(side = "top", fill = "x", pady = 4)

        Label(frame, text = text, font = self.title_font).pack(side = "left", padx = (0, 6))
        Frame(frame, height = 2, bd = 1, relief = "sunken").pack(side = "left", fill = "x", expand = True)

    def _build_slider(self, root: Frame, parameter: str, source_object: Orbit | Satellite, label: str, upper_lim: int, res: float = 1, lower_lim: int = 0) -> Scale:
        slider_name = f"_{parameter}_slider"
        self.__setattr__(
            slider_name,
            Scale(root, from_ = lower_lim, to = upper_lim, resolution = res, length = 195, orient = "horizontal",
                  command = partial(self._update_value, parameter, source_object), label = label, font = self.slider_font)
        )

        slider: Scale = self.__getattribute__(slider_name)
        init_value: float = round(getattr(source_object, parameter), 2) if parameter == "nu" else getattr(source_object, parameter)

        slider.set(init_value)
        slider.pack(side = "top", anchor = "nw")
        return slider

    def _update_value(self, parameter: str, source_object: Orbit | Satellite, new_val: str) -> None:
        new_val = float(new_val)
        setattr(source_object, parameter, new_val)

        self._orbit.update_orbital_properties()
        self._orbit.update_orbit_type()

        if parameter == "e":
            if new_val >= 1:
                t_asymp = round(self._orbit.t_asymp, 2)
                self._nu_slider.configure(from_ = -t_asymp, to = t_asymp)
                nu = self._sat.nu
                if nu < -t_asymp:
                    self._sat.nu = -t_asymp
                elif nu > t_asymp:
                    self._sat.nu = t_asymp

            else:
                self._nu_slider.configure(from_ = 0, to = 360)

        self._sat.update_satellite_properties()
        self._orbit_fig.redraw_orbit()

        for param_object, params in list(self._parameter_objects.items()):
            for param in params:
                self._update_display(param, param_object)


    def _build_display(self, frame: LabelFrame, parameter: str, source_object: Orbit | Satellite, display_str: str, units: str, row: int) -> None:
        var = StringVar(value = self._format_display_value(getattr(source_object, parameter), units))
        self.__setattr__(f"_{parameter}_str", var)

        name_label = Label(frame, text = display_str + ":", anchor = "w", font = self.slider_font)
        name_label.grid(row = row, column = 0, sticky = "w", padx = (0, 6))

        value_label = Label(frame, textvariable = var, anchor = "e", width = 13, font = self.slider_font)
        value_label.grid(row = row, column = 1, sticky = "ew", padx = (0, 6))

    def _update_display(self, parameter: str, source_object: Orbit | Satellite = None, value: float = None) -> None:
        new_value = value if value is not None else getattr(source_object, parameter)
        self.__getattribute__(
            f"_{parameter}_str"
        ).set(self._format_display_value(new_value, self.parameters[parameter][1]))

    def _format_display_value(self, value: float, units: str) -> str:
        if np.isclose(value, 0.00, rtol = 0.001):
            value = 0.00

        if np.isinf(value):
            return f"∞ {units}"

        elif np.isnan(value):
            return "n/a"

        elif units in ["km", "km²/s"]:
            return f"{value:6.0f} {units}"

        elif units in ["°", "km/s", "km²/s²"]:
            return f"{value:6.2f} {units}"