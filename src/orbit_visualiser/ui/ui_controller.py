from tkinter import Event
from typing import Callable
from orbit_visualiser.ui.input.input_panel_controller import InputController
from orbit_visualiser.ui.properties.properties_panel_controller import PropertiesController
from orbit_visualiser.ui.figure.orbit_figure_controller import OrbitFigureController
#from orbit_visualiser.ui.config.display_panel.display_panel_controller import DisplayController
from orbit_visualiser.ui.data_access import OrbitDataAccess
from orbit_visualiser.ui.ui_builder import UIBuilder


class UIController():


    def __init__(
            self,
            builder: UIBuilder,
            oda: OrbitDataAccess
    ):
        self._builder = builder

        self._figure_controller = OrbitFigureController(builder.figure_builder, oda)
        self._variables_controller = InputController(self._figure_controller, builder.input_builder, oda)
        self._properties_controller = PropertiesController(builder.properties_builder, oda)


    def validate_manual_input(
            self,
            variable: str,
            event: Event
    ) -> None:
        self._variables_controller.validate_manual_input(variable, event)
        self._properties_controller.update_display()

    def reset_state(self) -> Callable:
        return self._variables_controller.reset_state()

    def format_display_value(self, value: float | str, units: str | None) -> str:
        return self._properties_controller.format_display_value(value, units)

    def slider_changed(
            self,
            variable: str,
            input_type: str,
            new_val: str | float
    ) -> None:
        self._variables_controller.update_variable(
            variable,
            input_type,
            new_val
        )
        self._properties_controller.update_display()