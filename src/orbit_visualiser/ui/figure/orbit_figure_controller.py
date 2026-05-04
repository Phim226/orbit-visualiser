from mpl_toolkits.mplot3d.art3d import Line3D
from mpl_toolkits.mplot3d import Axes3D
from orbit_visualiser.ui.data_access import OrbitDataAccess
from orbit_visualiser.ui.figure.orbit_figure_builder import OrbitFigureBuilder

# TODO: Fix bug where scroll zoom doesn't register as changing the view so the native matplotlib home button has unexpected (and often undesirable) behaviour.
# TODO: Split into builder and controller
class OrbitFigureController():

    DISPLAY_TEXT_OFFSET = (1.5, 1.5)
    NUM_POINTS = 100_000

    def __init__(
            self,
            builder: OrbitFigureBuilder,
            oda: OrbitDataAccess
    ):
        self._builder = builder
        self._da = oda

    def redraw_orbit(self) -> None:
        x, y, z = self._da.get_orbit_data(OrbitFigureController.NUM_POINTS)
        self._builder.line.set_data_3d(x, y, z)

        self._builder.canvas.draw_idle()

    def redraw_satellite(self) -> None:
        x, y, z = self._da.get_sat_position()

        sat_point: Line3D = self._builder.satellite_point
        sat_point.set_xdata((x,))
        sat_point.set_ydata((y,))
        sat_point.set_3d_properties(zs = (z,))

        self._builder.canvas.draw_idle()

    def reset_axes(self) -> None:
        axis: Axes3D = self._builder.axis
        axis.set_xlim(-100_000, 100_000)
        axis.set_ylim(-100_000, 100_000)
        axis.set_zlim(-100_000, 100_000)

        self._builder.canvas.draw_idle()

    def plot_periapsis_point(self) -> None:
        axis: Axes3D = self._builder.axis
        self._rp_point, = axis.plot(
            self._da.satellite.orbit.radius_of_periapsis,
            0,
            ms = 3,
            marker = "o",
            zorder = 9,
            color = "#502BF2",
            label = "$r_p$"
        )
        self._rp_annotation = axis.annotate(
            "$r_p$",
            xy = (self._da.satellite.orbit.radius_of_periapsis, 0),
            xycoords = "data",
            xytext = OrbitFigureController.DISPLAY_TEXT_OFFSET,
            textcoords = "offset points"
        )