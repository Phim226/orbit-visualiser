from tkinter import Tk, Frame, Scale, Label, StringVar, LabelFrame, Button, Entry, DoubleVar, Checkbutton, IntVar
from tkinter.ttk import Separator
from typing import Callable
from functools import partial
import numpy as np
from orbit_visualiser.core import Orbit, Satellite, CentralBody
from orbit_visualiser.ui.config.display_panel.display_panel_builder import DisplayBuilder
from orbit_visualiser.ui.config.properties_panel.properties_panel_builder import PropertiesBuilder
from orbit_visualiser.ui.config.variables_panel.variables_panel_builder import VariablesBuilder
from orbit_visualiser.ui.common.builder import Builder
from orbit_visualiser.ui.common.specs import PropertySpec, VariableSpec

# TODO: Give option to show parameters on the plot (arrows/lines for vectors and distances etc).
# TODO: Split into variables, options and properties builders.
# TODO: Manage geometry of display options using rows/columns.

class OrbitConfigBuilder(Builder):

    def __init__(
            self, root: Tk,
            config_frame_placement: tuple[str],
            orbit: Orbit,
            central_body: CentralBody,
            satellite: Satellite
    ):
        self._root = root

        self._orbit = orbit
        self._central_body = central_body
        self._sat = satellite

        self._config_frame = Frame(root)
        self._config_frame.pack(
            side = config_frame_placement[0],
            anchor = config_frame_placement[1],
            padx = 8, pady = 6
        )

        self._options_frame = Frame(self._config_frame, padx = 2)
        self._options_frame.pack(side = "left", anchor = "n", pady = (2, 0))

        self._variables_builder = VariablesBuilder(self._options_frame, orbit, central_body, satellite)
        self._display_builder = DisplayBuilder(self._options_frame)
        self._properties_builder = PropertiesBuilder(self._config_frame, orbit, satellite)

    @property
    def variables_builder(self) -> VariablesBuilder:
        return self._variables_builder

    @property
    def display_builder(self) -> DisplayBuilder:
        return self._display_builder

    @property
    def properties_builder(self) -> PropertiesBuilder:
        return self._properties_builder

    def build(
            self,
            reset: Callable,
            validate_input: Callable,
            slider_changed: Callable,
            format_value: Callable
    ) -> None:
        self._variables_builder.build_variables_frame(reset, validate_input, slider_changed)
        self._display_builder.build_display_options_frame()

        sep = Separator(self._config_frame, orient = "vertical")
        sep.pack(side = "left", fill = "y", padx = 6, expand = True)

        self._properties_builder.build_properties_frame(format_value)