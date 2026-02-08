from tkinter import Event
from typing import Callable
from orbit_visualiser.core import Orbit, Satellite, CentralBody
from orbit_visualiser.ui.config.config_builder import OrbitConfigBuilder, OrbitConfigBuilder
from orbit_visualiser.ui.config.variables_panel.variables_panel_controller import VariablesController
from orbit_visualiser.ui.config.properties_panel.properties_panel_controller import PropertiesController
#from orbit_visualiser.ui.config.display_panel.display_panel_controller import DisplayController
from orbit_visualiser.ui.figure.orbit_figure import OrbitFigure


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