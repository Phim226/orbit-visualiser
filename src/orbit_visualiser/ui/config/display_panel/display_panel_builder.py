from tkinter import Tk, Frame, LabelFrame, Checkbutton, IntVar
from orbit_visualiser.ui.common.builder import Builder
from orbit_visualiser.core import Orbit, Satellite, CentralBody

# TODO: Give option to show parameters on the plot (arrows/lines for vectors and distances etc).
# TODO: Manage geometry of display options using rows/columns.

class DisplayBuilder(Builder):


    def __init__(self, options_frame: Frame):
        self._options_frame = options_frame

    def build_display_options_frame(self) -> None:
        display_options_frame: Frame = Frame(self._options_frame)
        self._display_options_frame = display_options_frame

        self._build_separator(display_options_frame, "Display options")

        # Build orbit options frame
        orbit_options_frame = LabelFrame(
            display_options_frame, bd = 1, relief = "sunken", text = "Orbit", font = self._subtitle_font
        )
        self._rp_display_var: IntVar = IntVar()
        rp_display_check: Checkbutton = Checkbutton(orbit_options_frame, text = "Periapsis", variable = self._rp_display_var)
        rp_display_check.pack(side = "top", anchor = "nw")
        orbit_options_frame.pack(side = "top", anchor = "nw", pady = (4, 0))

        # Build central body options frame
        central_body_options_frame = LabelFrame(
            display_options_frame, bd = 1, relief = "sunken", text = "Central body", font = self._subtitle_font
        )
        self._central_body_display_var: IntVar = IntVar(value = 1)
        central_body_display_check: Checkbutton = Checkbutton(central_body_options_frame, text = "Central body", variable = self._central_body_display_var)
        central_body_display_check.pack(side = "top", anchor = "nw")
        central_body_options_frame.pack(side = "top", anchor = "nw", pady = (4, 0))

        # Build satellite options frame
        satellite_options_frame = LabelFrame(
            display_options_frame, bd = 1, relief = "sunken", text = "Satellite", font = self._subtitle_font
        )
        self._radius_display_var: IntVar = IntVar()
        radius_display_check: Checkbutton = Checkbutton(satellite_options_frame, text = "Radius", variable = self._radius_display_var)
        radius_display_check.pack(side = "top", anchor = "nw")
        satellite_options_frame.pack(side = "top", anchor = "nw", pady = (4, 0))

        display_options_frame.pack(side = "top", anchor = "nw", pady = (2, 0), fill = "x", expand = True)