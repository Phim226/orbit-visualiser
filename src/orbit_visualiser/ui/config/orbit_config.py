from tkinter import Tk, Frame, Scale, Label, StringVar
from functools import partial
from orbit_visualiser.ui import OrbitFigure
from orbit_visualiser.core import Orbit

class OrbitConfigurer():

    orbital_parameters: dict[str, str] = {
        "a" : "Semi-major axis: ",
        "b" : "Semi-minor axis: ",
        "ra": "Radius of apoapsis: ",
        "p" : "Semi-parameter: "
    }

    def __init__(self, root: Tk, config_frame_placement: tuple[str], orbit_fig: OrbitFigure, orbit: Orbit):
        self._root = root

        self._orbit_fig = orbit_fig
        self._orbit = orbit

        self._config_frame = Frame(root)
        self._config_frame.pack(side = config_frame_placement[0], anchor = config_frame_placement[1])

    def build(self) -> None:
        self._build_eccentricity_slider()
        self._build_periapsis_slider()

        for parameter, display_str in list(self.orbital_parameters.items()):
            self._build_display(parameter, display_str)

    def _build_eccentricity_slider(self) -> None:
        self._e_slider = Scale(self._root, from_ = 0, to = 2, resolution = 0.01, length  = 150, orient = "horizontal",
                              command = partial(self._update_value, "e"), label = "Eccentricity")
        self._e_slider.set(self._orbit.e)
        self._e_slider.pack(side = "top", anchor = "nw")

    def _build_periapsis_slider(self) -> None:
        self._rp_slider = Scale(self._root, from_ = 0, to = 10, resolution = 0.01, length  = 150, orient = "horizontal",
                              command = partial(self._update_value, "rp"), label = "Radius of periapsis")
        self._rp_slider.set(self._orbit.rp)
        self._rp_slider.pack(side = "top", anchor = "nw")

    def _update_value(self, parameter: str, new_val: str) -> None:
        setattr(self._orbit, parameter, new_val)
        self._orbit_fig.redraw_orbit()

    def _build_display(self, parameter: str, display_str: str):
        self.__setattr__(
            f"_{parameter}_str",
            StringVar(value = f"{display_str}{getattr(self._orbit, parameter)}")
        )

        label = Label(self._root, textvariable = self.__getattribute__(f"_{parameter}_str"))
        label.pack(side = "top", anchor = "nw")

    #def _format_display(self, )