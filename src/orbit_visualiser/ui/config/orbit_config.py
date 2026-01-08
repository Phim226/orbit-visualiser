from tkinter import Tk, Frame, Scale
from orbit_visualiser.ui import OrbitFigure
from orbit_visualiser.core import Orbit

# TODO: Have options to display periapsis, apoapsis, semimajor/minor axes etc.
# TODO: Have option to adjust semimajor axis.
class OrbitConfigurer():


    def __init__(self, root: Tk, config_frame_placement: tuple[str], orbit_fig: OrbitFigure, orbit: Orbit):
        self._root = root

        self._orbit_fig = orbit_fig
        self._orbit = orbit

        self._config_frame = Frame(root)
        self._config_frame.pack(side = config_frame_placement[0], anchor = config_frame_placement[1])

    def build(self) -> None:
        self._build_eccentricity_slider()
        self._build_semimajor_axis_slider()

    def _build_eccentricity_slider(self) -> None:
        self._e_slider = Scale(self._root, from_ = 0, to = 2, resolution = 0.01, length  = 150, orient = "horizontal",
                              command = self._update_eccentricity, label = "Eccentricity")
        self._e_slider.set(self._orbit.e)
        self._e_slider.pack(side = "top", anchor = "nw")

    def _build_semimajor_axis_slider(self) -> None:
        self._a_slider = Scale(self._root, from_ = -10, to = 10, resolution = 0.01, length  = 150, orient = "horizontal",
                              command = self._update_semimajor_axis, label = "Semimajor axis")
        self._a_slider.set(self._orbit.a)
        self._a_slider.pack(side = "top", anchor = "nw")

    def _update_eccentricity(self, new_val: float) -> None:
        self._orbit.e = new_val
        self._a_slider.set(self._orbit.a)
        self._orbit_fig.redraw_orbit()

    def _update_semimajor_axis(self, new_val: float) -> None:
        self._orbit.a = new_val
        self._orbit_fig.redraw_orbit()