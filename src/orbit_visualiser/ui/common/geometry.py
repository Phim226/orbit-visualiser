from dataclasses import dataclass
from typing import Literal, LiteralString
from ttkbootstrap import Window

@dataclass
class DimensionsGeometry():
    width: int = 0
    height: int = 0

@dataclass
class PlaceableGeometry():
    x: int
    y: int

@dataclass
class FrameGeometry(DimensionsGeometry):
    side: Literal['left', 'right', 'top', 'bottom'] = "left"
    anchor: Literal['nw', 'n', 'ne', 'w', 'center', 'e', 'sw', 's', 'se'] = "nw"
    padx: int | tuple[int, int] = 0
    pady: int | tuple[int, int] = 0
    fill: Literal['none', 'x', 'y', 'both'] = "none"
    expand: bool = False

@dataclass
class ScrolledFrameGeometry(DimensionsGeometry):
    padding: int = 2

@dataclass
class EntryGeometry(DimensionsGeometry, PlaceableGeometry):
    pass

@dataclass
class SliderGeometry(PlaceableGeometry):
    pass

class GeometryManager():
    """
    Class for managing tkinter widget geometry and OS differences.
    """

    def __init__(self, os: LiteralString, root: Window):
        self._os = os
        self._os_win: bool = os.startswith("win")

        if self._os_win:
            root.geometry("1500x800")
            root.resizable(False, False)

    @property
    def parent_frame(self) -> tuple[FrameGeometry, FrameGeometry, FrameGeometry]:
        """
        The tuple of frame geometry of the 3 top level frames, input, figure and properties in that
        order.

        Returns
        -------
        tuple[FrameGeometry, FrameGeometry, FrameGeometry]
            Tuple of parent frame FrameGeometry
        """
        return (
            FrameGeometry(padx = 2, pady = (2, 0), fill = "y", expand = True),
            FrameGeometry(padx = 8, pady = 6),
            FrameGeometry(padx = 2, pady = (2, 0), fill = "y", expand = True),
        )

    @property
    def parent_scrollable(self) -> ScrolledFrameGeometry:
        """
        The geometry for the scrolled parent frames (input and properties)

        Returns
        -------
        ScrolledFrameGeometry
            Geometry for the parent ScrolledFrame
        """
        if self._os_win:
            return ScrolledFrameGeometry(width = 375, padding = 10)

        return ScrolledFrameGeometry(width = 315)

    @property
    def input_widgets(self) -> tuple[FrameGeometry, SliderGeometry, EntryGeometry]:
        """
        The geometry for the frame, slider and entry widgets for each individual input frame.

        Returns
        -------
        tuple[FrameGeometry, SliderGeometry, EntryGeometry]
            Geometry for the input frame
        """
        return (
            self._input_frame(),
            self._input_slider(),
            self._input_entry()
        )

    def _input_frame(self) -> FrameGeometry:
        if self._os_win:
            return FrameGeometry(side = "top", width = 375, height = 100, pady = 2)

        return FrameGeometry(side = "top", width = 290, height = 75, pady = 2)

    def _input_entry(self) -> EntryGeometry:
        if self._os_win:
            return EntryGeometry(x = 5, y = 30, width = 10)

        return EntryGeometry(x = 5, y = 20, width = 10)

    def _input_slider(self) -> SliderGeometry:
        if self._os_win:
            return SliderGeometry(x = 5, y = 70)

        return SliderGeometry(x = 5, y = 45)

    @property
    def figsize(self) -> tuple[int, int]:
        """
        The figure size of the orbit plot.

        Returns
        -------
        tuple[int, int]
            The tuple of ints representing the (x, y) figure size in inches
        """
        if self._os_win:
            return (5, 5)

        return (7, 7)