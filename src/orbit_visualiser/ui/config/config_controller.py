from tkinter import Entry, Event, messagebox, DoubleVar
from decimal import Decimal
import numpy as np
from orbit_visualiser.ui.common.specs import PropertySpec
from orbit_visualiser.ui.config.config_builder import OrbitConfigBuilder
from orbit_visualiser.ui.figure.orbit_figure import OrbitFigure
from orbit_visualiser.core import Orbit, Satellite, CentralBody

# TODO: Allow for temporary increase in slider scale when inputting manual values.
# TODO: Allow for fractional manual inputs.
# TODO: Remove any leading 0s from manual inputs.
# TODO: Refactor to a lazy/cached recalculation model. Currently everything is recalculated on every variable change.
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
        var_props = self._builder.variables_builder.variable_specs
        for name, value in list(var_props.items()):
            init_value = value.init_value
            getattr(self._builder.variables_builder, f"{name}_slider").set(init_value)

            entry: Entry = getattr(self._builder.variables_builder, f"{name}_entry")
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
        new_val = getattr(self._builder.variables_builder, f"{variable}_entry").get().strip()

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
                entry: Entry = getattr(self._builder.variables_builder, f"{variable}_entry")
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

    def slider_changed(
            self,
            variable: str,
            source_object: Orbit | Satellite | CentralBody,
            input_type: str,
            new_val: str | float
    ) -> None:
        self._variables_controller.update_variable(
            variable,
            source_object,
            input_type,
            new_val
        )
        self._properties_controller.update_display()