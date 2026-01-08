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
        slider_update = Scale(self._root, from_ = 0, to = 2, resolution = 0.01, orient = "horizontal",
                              command = self._orbit_fig.update_eccentricity, label = "Eccentricity")
        slider_update.set(self._orbit.e)
        slider_update.pack(side = "top")