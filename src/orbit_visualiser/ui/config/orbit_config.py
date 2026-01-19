from tkinter import Tk, Frame, Scale, Label, StringVar, LabelFrame, Button, Entry, Event, messagebox
from tkinter.ttk import Separator
from typing import Any
from functools import partial
from decimal import Decimal, ROUND_FLOOR
from math import copysign
import numpy as np
from orbit_visualiser.ui import OrbitFigure
from orbit_visualiser.core import Orbit, Satellite, CentralBody

# TODO: Give option to show parameters on the plot (arrows/lines for vectors and distances etc).
# TODO: Show the correct sign on the infinity symbol for x and y position.
# TODO: Allow for temporary increase in slider scale when inputting manual values.
# TODO: Refactor variable property dictionaries into dataclasses.
class OrbitConfigurer():

    title_font = ("Orbitron", 16, "bold")
    subtitle_font = ("Orbitron", 11, "normal")
    slider_font = ("Fira Mono", 9, "normal")

    orbital_parameters: dict[str, tuple[str]] = {
        "orbit_type" : ("Orbit type", None),
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
        "x" : ("x position", "km"),
        "y" : ("y position", "km"),
        "t" : ("Orbital period", "s"),
        "h" : ("Angular momentum", "km²/s"),
        "v" : ("Velocity", "km/s"),
        "v_azim" : ("Azimuthal velocity", "km/s"),
        "v_radial" : ("Radial velocity", "km/s"),
        "v_esc" : ("Escape velocity", "km/s"),
        "v_inf" : ("Excess velocity", "km/s"),
        "gam" : ("Flight angle", "°"),
        "eps" : ("Mechanical energy", "km²/s²"),
        "c3" : ("Characteristic energy", "km²/s²")
    }

    parameters: dict[str, tuple[str]] = orbital_parameters | satellite_parameters

    def __init__(self, root: Tk, config_frame_placement: tuple[str], orbit_fig: OrbitFigure, orbit: Orbit, central_body: CentralBody, satellite: Satellite):
        self._root = root

        self._orbit_fig = orbit_fig
        self._orbit = orbit
        self._central_body = central_body
        self._sat = satellite

        self._e_properties:  dict[str, Any] = {
            "name" : "Eccentricity",
            "object" : orbit,
            "init_value" : orbit.e,
            "slider_lims" : (0, 5),
            "decimal_places" : 3,
            "units" : None,
            "entry_pos" : (85, 4)
        }
        self._rp_properties: dict[str, Any] = {
            "name" : "Radius of periapsis",
            "object" : orbit,
            "init_value" : orbit.rp,
            "slider_lims" : (central_body.r + 1, 200_000),
            "decimal_places" : 0,
            "units" : "km",
            "entry_pos" : (160, 4)
        }
        self._mu_properties: dict[str, Any] = {
            "name" : "Gravitational parameter",
            "object" : central_body,
            "init_value" : central_body.mu,
            "slider_lims" : (1, 1_000_000),
            "decimal_places" : 0,
            "units" : "km³/s²",
            "entry_pos" : (198, 4)
        }
        self._nu_properties: dict[str, Any] = {
            "name" : "True anomaly",
            "object" : satellite,
            "init_value" : satellite.nu,
            "slider_lims" : (0, 360),
            "decimal_places" : 2,
            "units" : "°",
            "entry_pos" : (115, 4)
        }

        self._variable_properties: dict[str, dict] = {
            "e" : self._e_properties,
            "rp" : self._rp_properties,
            "mu" : self._mu_properties,
            "nu" : self._nu_properties
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

        # Build orbital geometry frame
        orbital_geom_frame = LabelFrame(var_frame, bd = 1, relief = "sunken", text = "Orbital geometry", font = self.subtitle_font)
        self._e_slider, self._e_entry = self._build_input_frame(orbital_geom_frame, "e", self._variable_properties["e"])
        self._rp_slider, self._rp_entry = self._build_input_frame(orbital_geom_frame, "rp", self._variable_properties["rp"])
        orbital_geom_frame.pack(side = "top", anchor = "nw", pady = (4, 0))

        # Building central body frame
        central_body_frame = LabelFrame(var_frame, bd = 1, relief = "sunken", text = "Central body", font = self.subtitle_font)
        self._mu_slider, self._mu_entry = self._build_input_frame(central_body_frame, "mu", self._variable_properties["mu"])
        central_body_frame.pack(side = "top", anchor = "nw", pady = (4, 0))

        # Build satellite frame
        sat_frame = LabelFrame(var_frame, bd = 1, relief = "sunken", text = "Satellite", font = self.subtitle_font)
        self._nu_slider, self._nu_entry = self._build_input_frame(sat_frame, "nu", self._variable_properties["nu"])
        sat_frame.pack(side = "top", anchor = "nw", pady = (4, 0))

        # Build reset button
        reset_button = Button(var_frame, text = "Reset", command = self._reset_state)
        reset_button.pack(side = "top", anchor = "nw", pady = (4, 0))

        var_frame.pack(side = "left", anchor = "n", pady = (2, 0))

    def _build_input_frame(self, root: Frame, parameter: str, param_props: dict[str, Any]) -> tuple[Scale, Entry]:
        obj = param_props["object"]
        units = param_props["units"]

        frame = Frame(root, width = 265, height = 60)

        slider = self._build_slider(
            frame,
            parameter,
            obj,
            f"{param_props["name"]}{"" if units is None else f" ({units})"} = ",
            param_props["slider_lims"],
            1/10**param_props["decimal_places"]
        )

        entry = Entry(frame, width = 10)
        entry.insert(0, f"{param_props["init_value"]: 0.{param_props["decimal_places"]}f}".strip())
        entry.bind("<Return>", partial(self._validate_manual_input, parameter, obj))
        x, y = param_props["entry_pos"]
        entry.place(x = x, y = y)

        frame.pack(side = "top", anchor = "nw", pady = 2)

        return slider, entry

    def _reset_state(self) -> None:
        for name, value in list(self._variable_properties.items()):
            init_value = value["init_value"]
            self.__getattribute__(f"_{name}_slider").set(init_value)

            entry = self.__getattribute__(f"_{name}_entry")
            entry.delete(0, 1000)
            entry.insert(0, f"{init_value: 0.{value["decimal_places"]}f}".strip())

            setattr(self._variable_properties[name]["object"], name, init_value)

        self._orbit.update_orbital_properties()
        self._orbit.update_orbit_type()
        self._sat.update_satellite_properties()

        self._orbit_fig.redraw_orbit()
        self._orbit_fig.redraw_satellite()
        self._orbit_fig.reset_axes()

    def _build_slider(self, root: Frame, parameter: str, source_object: Orbit | Satellite, label: str, lims: tuple[int], res: float) -> Scale:
        slider_name = f"_{parameter}_slider"
        self.__setattr__(
            slider_name,
            Scale(root, from_ = lims[0], to = lims[1], resolution = res, length = 260, orient = "horizontal",
                  command = partial(self._update_value, parameter, source_object, "slider"), label = label, font = self.slider_font)
        )

        slider: Scale = self.__getattribute__(slider_name)
        init_value: float = round(np.degrees(getattr(source_object, parameter)), 2) if parameter == "nu" else getattr(source_object, parameter)

        slider.set(init_value)
        slider.place(x = 0, y = 0, anchor = "nw")
        return slider

    def _validate_manual_input(self, parameter: str, source_object: Orbit | Satellite | CentralBody, event: Event) -> None:
        new_val = self.__getattribute__(f"_{parameter}_entry").get().strip()

        try:
            new_val_float = float(new_val)

            if new_val_float < 0 and parameter != "nu":
                raise ValueError

        except ValueError:
            messagebox.showwarning("Warning", "Invalid input")
            return

        # When e < 1 then the orbit is periodic, and so the true anomaly is as well.
        if parameter == "nu":
            if self._orbit.e < 1 and (new_val_float < 0 or new_val_float > 360):
                # float(new_val) will kill off any decimal points when new_val has extremely large absolute value. The Decimal class retains that information.
                # If the angle is negative then float(Decimal(new_val))%360 reduces it to (-360, 0), then + 360 to the range we want.
                new_val_float = (float(Decimal(new_val))%360 + 360)%360

            else:
                t_asymp = np.degrees(self._orbit.t_asymp)
                if new_val_float < -t_asymp:
                    new_val_float = -t_asymp
                elif new_val_float > t_asymp:
                    new_val_float = t_asymp

        self._update_value(parameter, source_object, "entry", new_val_float)

    def _update_value(self, parameter: str, source_object: Orbit | Satellite, input_type: str, new_val: str | float) -> None:
        new_val = float(new_val)

        # This if-elif block lets the sliders and manual inputs update one another.
        if input_type == "slider":
            entry: Entry = self.__getattribute__(f"_{parameter}_entry")
            entry.delete(0, 1000)
            entry.insert(0, f"{new_val: 0.{self._variable_properties[parameter]["decimal_places"]}f}".strip())
        elif input_type == "entry":
            self.__getattribute__(f"_{parameter}_slider").set(new_val)

        if parameter == "nu":
            new_val = np.deg2rad(float(new_val))

        setattr(source_object, parameter, new_val)


        self._orbit.update_orbital_properties()
        self._orbit.update_orbit_type()

        # The value of the eccentricity determines the range of possible true anomaly values, which this if block checks for.
        if parameter == "e":
            if new_val >= 1:
                t_asymp = self._orbit.t_asymp
                t_asymp_slider_lim = round(np.degrees(t_asymp), 2)
                self._nu_slider.configure(from_ = -t_asymp_slider_lim, to = t_asymp_slider_lim)
                nu = self._sat.nu
                if nu < -t_asymp:
                    self._sat.nu = -t_asymp
                elif nu > t_asymp:
                    self._sat.nu = t_asymp

            else:
                self._nu_slider.configure(from_ = 0, to = 360)

        self._sat.update_satellite_properties()

        self._orbit_fig.redraw_orbit()
        self._orbit_fig.redraw_satellite()

        for param_object, params in list(self._parameter_objects.items()):
            for param in params:
                self._update_display(param, param_object)

    def _update_display(self, parameter: str, source_object: Orbit | Satellite = None, value: float = None) -> None:
        new_value = value if value is not None else getattr(source_object, parameter)
        if self.parameters[parameter][1] == "°":
            new_value = np.degrees(new_value)

        self.__getattribute__(
            f"_{parameter}_str"
        ).set(self._format_display_value(new_value, self.parameters[parameter][1]))

    def _format_display_value(self, value: float | str, units: str | None) -> str:
        if units is None:
            return value.capitalize()

        if np.isclose(value, 0.00, rtol = 0.001):
            value = 0.00

        if np.isneginf(value):
            return f"-∞ {units}"

        elif np.isinf(value):
            return f"∞ {units}"

        elif np.isnan(value):
            return "n/a"

        elif units in ["km", "km²/s", "s"]:
            return f"{value:6.0f} {units}"

        elif units in ["°", "km/s", "km²/s²"]:
            return f"{value:6.2f} {units}"

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

    def _build_display(self, frame: LabelFrame, parameter: str, source_object: Orbit | Satellite, display_str: str, units: str, row: int) -> None:
        var = StringVar(value = self._format_display_value(getattr(source_object, parameter), units))
        self.__setattr__(f"_{parameter}_str", var)

        name_label = Label(frame, text = display_str + ":", anchor = "w", font = self.slider_font)
        name_label.grid(row = row, column = 0, sticky = "w", padx = (0, 6))

        value_label = Label(frame, textvariable = var, anchor = "e", width = 13, font = self.slider_font)
        value_label.grid(row = row, column = 1, sticky = "ew", padx = (0, 6))

    def _build_separator(self, root: Frame, text: str) -> None:
        frame = Frame(root)
        frame.pack(side = "top", fill = "x", pady = 4)

        Label(frame, text = text, font = self.title_font).pack(side = "left", padx = (0, 6))
        Frame(frame, height = 2, bd = 1, relief = "sunken").pack(side = "left", fill = "x", expand = True)