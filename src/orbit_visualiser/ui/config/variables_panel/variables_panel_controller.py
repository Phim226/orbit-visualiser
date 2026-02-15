from tkinter import Entry, Event, messagebox, DoubleVar
from decimal import Decimal
import numpy as np
from orbit_visualiser.ui.figure.orbit_figure import OrbitFigure
from orbit_visualiser.ui.config.variables_panel.variables_panel_builder import  VariablesBuilder
from orbit_visualiser.core import Orbit, Satellite, CentralBody
from orbit_visualiser.core.satellite import NewSatellite
from orbit_visualiser.core.neworbit import NewOrbit
from orbit_visualiser.core.astrodynamics.keplerian.elements import asymptote_anomaly
from orbit_visualiser.core.astrodynamics.types import OrbitType

# TODO: Allow for temporary increase in slider scale when inputting manual values.
# TODO: Allow for fractional manual inputs.
# TODO: Remove any leading 0s from manual inputs.
# TODO: Refactor to a lazy/cached recalculation model. Currently everything is recalculated on every variable change.

class VariablesController():


    def __init__(
            self,
            figure: OrbitFigure,
            builder: VariablesBuilder,
            satellite: NewSatellite
    ):
        self._orbit_fig = figure
        self._builder = builder

        self._satellite = satellite

    def reset_state(self) -> None:
        init_values = []
        var_props = self._builder.variable_specs
        for name, value in list(var_props.items()):
            init_value = value.init_value
            getattr(self._builder, f"{name}_slider").set(init_value)

            entry: Entry = getattr(self._builder, f"{name}_entry")
            entry.delete(0, 1000)
            entry.insert(0, f"{init_value: 0.{value.decimal_places}f}".strip())

            init_values.append(init_value)

        self._update_satellite_state(*init_values)

        self._orbit_fig.redraw_orbit()
        self._orbit_fig.redraw_satellite()
        self._orbit_fig.reset_axes()

    def validate_manual_input(
            self,
            variable: str,
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
            if self._satellite.orbit.eccentricity < 1 and (new_val_float < 0 or new_val_float > 360):
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
                t_asymp = np.degrees(self._satellite.orbit.asymptote_anomaly)
                if new_val_float < -t_asymp:
                    new_val_float = -t_asymp
                elif new_val_float > t_asymp:
                    new_val_float = t_asymp

        self.update_variable(variable, "entry", new_val_float)

    def update_variable(
            self,
            variable: str,
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

        new_values = {
            "e": self._builder.e_var.get(),
            "rp": self._builder.rp_var.get(),
            "mu": self._builder.mu_var.get(),
            "nu": np.deg2rad(self._builder.nu_var.get()),
        }
        new_values[variable] = new_val

        # The value of the eccentricity determines the range of possible true anomaly values, which
        # this if block checks for.
        t_asymp = np.nan
        if variable == "e":
            if new_val >= 1:
                t_asymp = asymptote_anomaly(OrbitType.HYPERBOLIC, new_val)
                t_asymp_slider_lim = round(np.degrees(t_asymp), 2)
                self._builder.nu_slider.configure(from_ = -t_asymp_slider_lim, to = t_asymp_slider_lim)

                nu = self._satellite.orbit.true_anomaly
                if nu < -t_asymp:
                    new_values["nu"] = -t_asymp
                elif nu > t_asymp:
                    new_values["nu"] = t_asymp
            else:
                self._builder.nu_slider.configure(from_ = 0, to = 360)

        self._update_satellite_state(*new_values.values(), t_asymp)

        self._orbit_fig.redraw_orbit()
        self._orbit_fig.redraw_satellite()

    def _update_satellite_state(self, e: float, rp: float, mu: float, nu: float, asymptote_anomaly: float = np.nan) -> None:
        orbit = NewOrbit.from_orbital_elements(e, rp, mu, nu, asymptote_anomaly)
        self._satellite.position = orbit.position
        self._satellite.velocity = orbit.velocity
        self._satellite.central_body.mu = mu