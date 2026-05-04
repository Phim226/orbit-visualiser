from tkinter import Frame
from ttkbootstrap.scrolled import ScrolledFrame
from typing import Callable
from orbit_visualiser.ui.figure.orbit_figure_builder import OrbitFigureBuilder
from orbit_visualiser.ui.properties.properties_panel_builder import PropertiesBuilder
from orbit_visualiser.ui.input.input_panel_builder import InputBuilder
from orbit_visualiser.ui.data_access import OrbitDataAccess

class UIBuilder():
    INPUT_GEOMETRY: dict[str, str | int | bool] = {
        "side": "left",
        "anchor": "nw",
        "padx": 2,
        "pady": (2, 0),
        "fill": "y",
        "expand": True
    }
    FIGURE_GEOMETRY: dict[str, str | int | bool] = {
        "side": "left",
        "anchor": "nw",
        "padx": 8,
        "pady": 6,
        "fill": "both",
        "expand": True
    }
    PROPS_GEOMETRY: dict[str, str | int | bool] = {
        "side": "left",
        "anchor": "nw",
        "padx": 0,
        "pady": (2, 0),
        "fill": "y",
        "expand": True
    }

    def __init__(self, root: Frame, oda: OrbitDataAccess):
        self._root = root
        self._oda = oda

        self._input_frame = self._build_frame(UIBuilder.INPUT_GEOMETRY)
        self._figure_frame = self._build_frame(UIBuilder.FIGURE_GEOMETRY, scrollable = False)
        self._properties_frame = self._build_frame(UIBuilder.PROPS_GEOMETRY)

        self._input_builder = InputBuilder(self._input_frame, oda)
        self._figure_builder = OrbitFigureBuilder(self._figure_frame, oda)
        self._properties_builder = PropertiesBuilder(self._properties_frame, oda)

    @property
    def input_builder(self) -> InputBuilder:
        return self._input_builder

    @property
    def figure_builder(self) -> OrbitFigureBuilder:
        return self._figure_builder

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
        self._input_builder.build_input_frame(reset, validate_input, slider_changed)
        self._figure_builder.build()
        self._properties_builder.build_properties_frame(format_value)

    def _build_frame(self, frame_geom: dict[str, str | int | bool], scrollable: bool = True) -> Frame | ScrolledFrame:
        frame = (ScrolledFrame(self._root, width = 315, autohide = True)
                if scrollable else Frame(self._root, bd = 1))
        frame.pack(side = frame_geom["side"], anchor = frame_geom["anchor"], padx = frame_geom["padx"],
                   pady = frame_geom["pady"], fill = frame_geom["fill"], expand = frame_geom["expand"]
                   )

        return frame


