from tkinter import Entry, Event, messagebox, DoubleVar
from decimal import Decimal
import numpy as np
from orbit_visualiser.ui import OrbitFigure, OrbitConfigBuilder, ParameterSpec
from orbit_visualiser.core import Orbit, Satellite, CentralBody

# TODO: Show the correct sign on the infinity symbol for x and y position.
# TODO: Allow for temporary increase in slider scale when inputting manual values.
# TODO: Allow for fractional manual inputs.
# TODO: Remove any leading 0s from manual inputs.
# TODO: Refactor to a lazy/cached recalculation model. Currently everything is recalculated on every variable change.
# TODO: Show correct sign on infinity symbol for mean anomaly and time since periapsis for open trajectories.
# TODO: Properly display very large values in the properties panel without them being cut off.

class OrbitConfigController():


    def __init__(
            self,
            figure: OrbitFigure,
            builder: OrbitConfigBuilder,
            orbit: Orbit,
            satellite: Satellite,
            central_body: CentralBody
    ):
        self._orbit_fig: OrbitFigure = figure
        self._builder: OrbitConfigBuilder = builder

        self._orbit: Orbit = orbit
        self._sat: Satellite = satellite
        self._central_body: CentralBody = central_body


    def reset_state(self) -> None:
        var_props = self._builder.variable_specs
        for name, value in list(var_props.items()):
            init_value = value.init_value
            getattr(self._builder, f"{name}_slider").set(init_value)

            entry: Entry = getattr(self._builder, f"{name}_entry")
            entry.delete(0, 1000)
            entry.insert(0, f"{init_value: 0.{value.decimal_places}f}".strip())

            setattr(var_props[name].obj, name, init_value)

        self._orbit.update_orbital_properties()
        self._sat.update_satellite_properties()

        self._orbit_fig.redraw_orbit()
        self._orbit_fig.redraw_satellite()
        self._orbit_fig.reset_axes()

    def validate_manual_input(
            self,
            variable: str,
            source_object: Orbit | Satellite | CentralBody,
            event: Event
    ) -> None:
        new_val = getattr(self._builder, f"{variable}_entry").get().strip()

        try:
            new_val_float = float(new_val)

            if new_val_float < 0 and variable != "nu":
                raise ValueError

        except ValueError:
            messagebox.showwarning("Warning", "Invalid input")
            return

        # When e < 1 then the orbit is periodic, and so the true anomaly is as well.
        if variable == "nu":
            if self._orbit.e < 1 and (new_val_float < 0 or new_val_float > 360):
                # float(new_val) will kill off any decimal points when new_val has extremely large
                # absolute value (around 16 digits due to limitations of 64bit double precision
                # for python floats). The Decimal class retains that information. If the angle is
                # negative then Decimal(new_val)%360 reduces it to (-360, 0), then + 360 to the range we want.
                new_val_float = (Decimal(new_val)%360 + 360)%360
                entry: Entry = getattr(self._builder, f"{variable}_entry")
                entry.delete(0, 1000)
                entry.insert(
                    0,
                    f"{new_val_float: 0.{self._builder.variable_specs[variable].decimal_places}f}".strip()
                )

            else:
                t_asymp = np.degrees(self._orbit.t_asymp)
                if new_val_float < -t_asymp:
                    new_val_float = -t_asymp
                elif new_val_float > t_asymp:
                    new_val_float = t_asymp

        self.update_value(variable, source_object, "entry", new_val_float)

    def update_value(
            self,
            variable: str,
            source_object: Orbit | Satellite | CentralBody,
            input_type: str,
            new_val: str | float
    ) -> None:
        new_val = float(new_val)

        # This if-elif block lets the sliders and manual inputs update one another.
        if input_type == "slider":
            entry: Entry = getattr(self._builder, f"{variable}_entry")
            entry.delete(0, 1000)
            entry.insert(
                0,
                f"{new_val: 0.{self._builder.variable_specs[variable].decimal_places}f}".strip()
            )

        elif input_type == "entry":
            slider_var: DoubleVar = getattr(self._builder, f"{variable}_var")
            slider_var.set(new_val)

        if variable == "nu":
            new_val = np.deg2rad(float(new_val))

        setattr(source_object, variable, new_val)


        self._orbit.update_orbital_properties()

        # The value of the eccentricity determines the range of possible true anomaly values, which
        # this if block checks for.
        if variable == "e":
            if new_val >= 1:
                t_asymp = self._orbit.t_asymp
                t_asymp_slider_lim = round(np.degrees(t_asymp), 2)
                self._builder.nu_slider.configure(from_ = -t_asymp_slider_lim, to = t_asymp_slider_lim)
                nu = self._sat.nu
                if nu < -t_asymp:
                    self._sat.nu = -t_asymp
                elif nu > t_asymp:
                    self._sat.nu = t_asymp

            else:
                self._builder.nu_slider.configure(from_ = 0, to = 360)

        self._sat.update_satellite_properties()

        self._orbit_fig.redraw_orbit()
        self._orbit_fig.redraw_satellite()

        for property_object, properties in list(self._builder.property_specs_by_object.items()):
            for property in properties:
                self._update_display(property, property_object)

    def _update_display(
            self,
            property: str,
            source_object: Orbit | Satellite = None,
            value: float = None
    ) -> None:
        spec: PropertySpec = self._builder.properties[property]
        new_value = value if value is not None else spec.getter(source_object)
        unit = spec.units
        if unit is not None and "°" in unit:
            new_value = np.degrees(new_value)

        getattr(self._builder, f"{property}_str").set(self.format_display_value(new_value, unit))

    def format_display_value(self, value: float | str, units: str | None) -> str:
        if units is None:
            return value.capitalize()

        if np.isclose(value, 0):
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

        elif units in ["°/s"]:
            return f"{value:6.6f} {units}"