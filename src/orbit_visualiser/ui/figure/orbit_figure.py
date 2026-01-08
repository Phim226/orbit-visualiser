from tkinter import Tk, Frame
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from orbit_visualiser.core import Orbit

# TODO: Figure out nice way of displaying parabolic orbits (maybe have option to make some parameters constant).
# TODO: Display point at pericenter (or the origin of the grid in this case since we are in the perifocal frame).
class OrbitFigure():

    def __init__(self, root: Tk, figure_frame_placement: tuple[str], orbit: Orbit):
        self._root = root

        self._orbit = orbit

        self._figure_frame: Frame = Frame(root)
        self._figure_frame.pack(side = figure_frame_placement[0], anchor = figure_frame_placement[1])

    def build(self) -> None:
        self._create_figure()
        self._configure_axes()
        self._initialise_plot()

        self._build_canvas()
        self._build_toolbar()

    def _create_figure(self) -> None:
        self._fig = Figure(figsize = (5, 4), dpi = 100)
        self._fig.subplots_adjust(left = 0, right = 1.0, bottom = 0, top = 1.0)
        self._ax = self._fig.add_subplot()

    def _configure_axes(self) -> None:
        self._ax.spines['left'].set_position(('data', 0))
        self._ax.spines['bottom'].set_position(('data', 0))

        self._ax.spines['right'].set_color('none')
        self._ax.spines['top'].set_color('none')

        self._ax.xaxis.set_ticks_position('bottom')
        self._ax.yaxis.set_ticks_position('left')
        self._ax.set_xlim(-10, 10)
        self._ax.set_ylim(-10, 10)

    def _initialise_plot(self) -> None:
        t = self._orbit.orbital_angles()
        orbit_eq = self._orbit.orbit_eq
        self._line, = self._ax.plot(orbit_eq.x(t) , orbit_eq.y(t))

    def _build_canvas(self) -> None:
        self._canvas = FigureCanvasTkAgg(self._fig, master = self._figure_frame)
        self._canvas.draw()
        self._canvas.get_tk_widget().pack(side = "top", fill = "both", expand = True)

    def _build_toolbar(self) -> None:
        toolbar = NavigationToolbar2Tk(self._canvas, self._figure_frame, pack_toolbar = False)
        toolbar.update()
        toolbar.pack(side = "bottom", fill = "x")

    def redraw_orbit(self) -> None:
        t = self._orbit.orbital_angles()
        orbit_eq = self._orbit.orbit_eq
        self._line.set_data(orbit_eq.x(t), orbit_eq.y(t))

        self._canvas.draw()

