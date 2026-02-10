from tkinter import Tk, Frame
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.patches import Circle
import numpy as np
from numpy.typing import NDArray
from orbit_visualiser.core import Orbit, CentralBody, Satellite

# TODO: Fix bug where scroll zoom doesn't register as changing the view so the native matplotlib home button has unexpected (and often undesirable) behaviour.
class OrbitFigure():

    DISPLAY_TEXT_OFFSET = (1.5, 1.5)
    NUM_POINTS = 100_000

    def __init__(
            self,
            root: Tk,
            figure_frame_placement: tuple[str],
            orbit: Orbit,
            central_body: CentralBody,
            satellite : Satellite
    ):
        self._root = root

        self._orbit = orbit
        self._body = central_body
        self._sat = satellite

        self._figure_frame: Frame = Frame(root)
        self._figure_frame.pack(
            side = figure_frame_placement[0], anchor = figure_frame_placement[1],
            padx = 8, pady = 6, fill = "both", expand = True
        )

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
        axis_colour = "#4D4D4DFF"

        self._ax.spines['left'].set_position(('data', 0))
        self._ax.spines['bottom'].set_position(('data', 0))

        self._ax.spines['left'].set_color(axis_colour)
        self._ax.spines['bottom'].set_color(axis_colour)
        self._ax.spines['right'].set_color('none')
        self._ax.spines['top'].set_color('none')

        self._ax.xaxis.set_ticks_position('bottom')
        self._ax.yaxis.set_ticks_position('left')
        self._ax.tick_params(colors = axis_colour)

        self._ax.set_xlim(-100_000, 100_000)
        self._ax.set_ylim(-100_000, 100_000)

        self._ax.text(
            0.98, 0.02,
            r"$\mathrm{km}$",
            transform = self._ax.transAxes,
            ha = "right",
            va = "bottom",
            fontsize = 9,
            color = "gray"
        )

        self._ax.tick_params(labelsize = 8)

    def _initialise_plot(self) -> None:
        # Plot the initial orbit
        nu = self._get_anomaly_data(self._orbit.orbital_angles())
        x, y = self._sat.pos_pf_eq(nu)
        self._line, = self._ax.plot(x, y, color = "#2F2F2F", alpha = 0.5, linewidth = 1.5)

        # Plot the central body
        self._ax.add_patch(
            Circle((0, 0), radius = self._body.r, fill = True, zorder = 10,
                   facecolor = "#4C6A92", edgecolor = "#3C5474")
        )

        # Plot the satellite
        self._sat_point, = self._ax.plot(
            self._orbit.rp, 0, ms = 5, marker = "o", zorder = 10, color = "#F28E2B"
        )

        self.plot_periapsis_point()

        f = self._zoom_factory(self._ax, 1.1)

    def _build_canvas(self) -> None:
        self._canvas = FigureCanvasTkAgg(self._fig, master = self._figure_frame)
        self._canvas.draw()
        self._canvas.get_tk_widget().pack(side = "top", fill = "both", expand = True)

    def _build_toolbar(self) -> None:
        toolbar = NavigationToolbar2Tk(self._canvas, self._figure_frame, pack_toolbar = False)
        toolbar.update()
        toolbar.pack(side = "bottom", fill = "x")

    def _get_anomaly_data(self, angles: tuple[float]) -> NDArray[np.float64]:
        lower_lim, upper_lim = angles

        # Check if angles represent closed or open orbit (open orbits will have negative angles).
        # If the orbit is open then we need a small offset so the plot doesn't evaluate to infinity
        # and cause runtime errors or unusual graphical artifacts.
        delta = 0.0001 if lower_lim < 0 else 0
        return np.linspace(lower_lim + delta, upper_lim - delta, OrbitFigure.NUM_POINTS)

    def redraw_orbit(self) -> None:
        nu = self._get_anomaly_data(self._orbit.orbital_angles())
        x, y = self._sat.pos_pf_eq(nu)
        self._line.set_data(x, y)

        self._canvas.draw_idle()

    def redraw_satellite(self) -> None:
        x, y = self._sat.pos_pf
        self._sat_point.set_data((x,), (y,))

        self._canvas.draw_idle()

    def reset_axes(self) -> None:
        self._ax.spines['left'].set_position(('data', 0))
        self._ax.spines['bottom'].set_position(('data', 0))
        self._ax.set_xlim(-100_000, 100_000)
        self._ax.set_ylim(-100_000, 100_000)

        self._canvas.draw_idle()

    def plot_periapsis_point(self) -> None:

        self._rp_point, = self._ax.plot(
            self._orbit.rp, 0, ms = 3, marker = "o", zorder = 9, color = "#502BF2", label = "$r_p$"
        )
        self._rp_annotation = self._ax.annotate(
            "$r_p$",
            xy = (self._orbit.rp, 0),
            xycoords = "data",
            xytext = OrbitFigure.DISPLAY_TEXT_OFFSET,
            textcoords = "offset points"
        )

    @staticmethod
    def _zoom_factory(ax: Axes, base_scale = 2.):
        def zoom_fun(event):
            # get the current x and y limits
            cur_xlim = ax.get_xlim()
            cur_ylim = ax.get_ylim()
            # set the range
            cur_xrange = (cur_xlim[1] - cur_xlim[0])*.5
            cur_yrange = (cur_ylim[1] - cur_ylim[0])*.5
            plot_centre_x = (cur_xlim[1] + cur_xlim[0])/2
            plot_centre_y = (cur_ylim[1] + cur_ylim[0])/2
            if event.button == 'up':
                # deal with zoom in
                scale_factor = 1/base_scale
            elif event.button == 'down':
                # deal with zoom out
                scale_factor = base_scale
            else:
                # deal with something that should never happen
                scale_factor = 1
                print(event.button)
            # set new limits
            ax.set_xlim([plot_centre_x - cur_xrange*scale_factor,
                        plot_centre_x + cur_xrange*scale_factor])
            ax.set_ylim([plot_centre_y - cur_yrange*scale_factor,
                        plot_centre_y + cur_yrange*scale_factor])
            ax.figure.canvas.draw_idle() # force re-draw the next time the GUI refreshes

        fig = ax.get_figure() # get the figure of interest
        # attach the call back
        fig.canvas.mpl_connect('scroll_event',zoom_fun)

        #return the function
        return zoom_fun

