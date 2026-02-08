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
        self._orbit_fig = figure
        self._builder = builder

        self._orbit = orbit
        self._sat = satellite

        self._variables_controller = VariablesController(figure, builder.variables_builder, orbit, satellite, central_body)
        self._properties_controller = PropertiesController(builder.properties_builder)


    def validate_manual_input(
            self,
            variable: str,
            source_object: Orbit | Satellite | CentralBody,
            event: Event
    ) -> None:
        return self._variables_controller.validate_manual_input(variable, source_object, event)

    def reset_state(self) -> Callable:
        return self._variables_controller.reset_state()

    def format_display_value(self, value: float | str, units: str | None) -> str:
        return self._properties_controller.format_display_value(value, units)

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