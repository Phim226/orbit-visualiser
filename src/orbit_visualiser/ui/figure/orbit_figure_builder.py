from ttkbootstrap import Frame
import matplotlib as mpl
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases import NavigationToolbar2
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3D
import numpy as np
from orbit_visualiser.ui.data_access import OrbitDataAccess

# TODO: Fix bug where scroll zoom doesn't register as changing the view so the native matplotlib home button has unexpected (and often undesirable) behaviour.
# TODO: Split into builder and controller
class OrbitFigureBuilder():

    DISPLAY_TEXT_OFFSET = (1.5, 1.5)
    NUM_POINTS = 100_000

    def __init__(
            self,
            figure_frame: Frame,
            oda : OrbitDataAccess
    ):
        self._figure_frame = figure_frame
        self._da = oda

    @property
    def line(self) -> Line3D:
        return self._line

    @property
    def satellite_point(self) -> Line3D:
        return self._sat_point

    @property
    def canvas(self) -> FigureCanvasTkAgg:
        return self._canvas

    @property
    def axis(self) -> Axes3D:
        return self._ax

    def build(self) -> None:
        self._create_figure()
        self._configure_axes()
        self._configure_figure_parameters()
        self._initialise_plot()

        self._build_canvas()
        self._build_toolbar()

    def _create_figure(self) -> None:
        self._fig = Figure(figsize = (7, 6), dpi = 100)
        self._fig.subplots_adjust(left = 0, right = 1.1, bottom = -0.1, top = 1.1)
        self._ax: Axes3D = self._fig.add_subplot(projection = "3d")
        self._ax.set_aspect("equal", adjustable = "datalim")

    def _configure_axes(self) -> None:
        axis_colour = "#4D4D4DFF"

        self._ax.xaxis.set_ticks_position('lower')
        self._ax.yaxis.set_ticks_position('lower')
        self._ax.zaxis.set_ticks_position('lower')
        self._ax.tick_params(colors = axis_colour)

        self._ax.set_xlim(-100_000, 100_000)
        self._ax.set_ylim(-100_000, 100_000)
        self._ax.set_zlim(-100_000, 100_000)

    @staticmethod
    def _configure_figure_parameters() -> None:
        mpl.rcParams['axes3d.mouserotationstyle'] = 'azel'

        NavigationToolbar2.toolitems = (
            ('Home', 'Reset original view', 'home', 'home'),
            ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'),
            ('Save', 'Save the figure', 'filesave', 'save_figure'),
        )

    def _plot_central_body(self) -> None:
        u, v = np.meshgrid(np.linspace(0, 2*np.pi, 25), np.linspace(0, np.pi, 25))
        r = self._da.satellite.central_body.r
        x = r*np.cos(u)*np.sin(v)
        y = r*np.sin(u)*np.sin(v)
        z = r*np.cos(v)
        self._ax.plot_wireframe(x, y, z, zorder = 10, edgecolor = "#3C5474")

    def _initialise_plot(self) -> None:
        # Plot the initial orbit
        x, y, z = self._da.get_orbit_data(OrbitFigureBuilder.NUM_POINTS)
        self._line: Line3D = self._ax.plot(x, y, z, color = "#2F2F2F", alpha = 0.5, linewidth = 1.5)[0]

        # Plot the central body
        self._plot_central_body()

        # Plot the satellite
        self._sat_point: Line3D = self._ax.plot(
            self._da.satellite.orbit.radius_of_periapsis, 0, 0, ms = 5, marker = "o", zorder = 10, color = "#F28E2B"
        )[0]

        #self.plot_periapsis_point()

        self._zoom_factory(self._ax, 1.1)

    def _build_canvas(self) -> None:
        self._canvas = FigureCanvasTkAgg(self._fig, master = self._figure_frame)
        self._canvas.draw()
        self._canvas.get_tk_widget().pack(side = "top", fill = "both", expand = True)

    def _build_toolbar(self) -> None:
        toolbar = NavigationToolbar2Tk(self._canvas, self._figure_frame, pack_toolbar = False)
        toolbar.update()
        toolbar.pack(side = "bottom", fill = "x")

    @staticmethod
    def _zoom_factory(ax: Axes3D, base_scale = 2.):
        def zoom_fun(event):
            cur_xlim = ax.get_xlim()
            cur_ylim = ax.get_ylim()
            cur_zlim = ax.get_zlim()

            cur_xrange = (cur_xlim[1] - cur_xlim[0])*.5
            cur_yrange = (cur_ylim[1] - cur_ylim[0])*.5
            cur_zrange = (cur_zlim[1] - cur_zlim[0])*.5
            plot_centre_x = (cur_xlim[1] + cur_xlim[0])/2
            plot_centre_y = (cur_ylim[1] + cur_ylim[0])/2
            plot_centre_z = (cur_zlim[1] + cur_zlim[0])/2

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

            ax.set_xlim([plot_centre_x - cur_xrange*scale_factor,
                        plot_centre_x + cur_xrange*scale_factor])
            ax.set_ylim([plot_centre_y - cur_yrange*scale_factor,
                        plot_centre_y + cur_yrange*scale_factor])
            ax.set_zlim([plot_centre_z - cur_zrange*scale_factor,
                        plot_centre_z + cur_zrange*scale_factor])
            ax.figure.canvas.draw_idle()

        fig = ax.get_figure()
        fig.canvas.mpl_connect('scroll_event',zoom_fun)

        return zoom_fun

