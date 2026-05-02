import sys
from tkinter import Tk
from orbit_visualiser.core import Orbit, Satellite, CentralBody
from orbit_visualiser.ui import OrbitFigure, OrbitConfigBuilder, OrbitConfigController
from orbit_visualiser.ui.common.presets import initial_config
from orbit_visualiser.ui.data_access import OrbitDataAccess

# TODO: Redesign layout geometry

class OrbitVisualiser():

    FIGURE_GEOMETRY = ("left", "nw")
    CONFIG_GEOMETRY = ("right", "ne")

    def __init__(self, root: Tk):
        root.title("2D Orbit Visualiser")

        if sys.platform.startswith("win"):
            root.state("zoomed")
        else:
            root.state("normal")

        da: OrbitDataAccess = self._initialise_orbit_objects()

        orbit_figure: OrbitFigure = OrbitFigure(root, OrbitVisualiser.FIGURE_GEOMETRY, da)
        orbit_figure.build()

        orbit_builder: OrbitConfigBuilder = OrbitConfigBuilder(root, OrbitVisualiser.CONFIG_GEOMETRY, da)
        orbit_controller: OrbitConfigController = OrbitConfigController(orbit_figure, orbit_builder, da)

        orbit_builder.build(
            orbit_controller.reset_state,
            orbit_controller.validate_manual_input,
            orbit_controller.slider_changed,
            orbit_controller.format_display_value
        )

    def _initialise_orbit_objects(self) -> OrbitDataAccess:
        orbit: Orbit = Orbit.from_orbital_elements(
            initial_config.eccentricity,
            initial_config.radius_of_periapsis,
            initial_config.true_anomaly,
            initial_config.right_ascension_of_the_ascending_node,
            initial_config.inclination,
            initial_config.argument_of_periapsis,
            initial_config.gravitational_parameter
        )

        central_body: CentralBody = CentralBody(
            initial_config.gravitational_parameter,
            initial_config.radius
        )
        satellite: Satellite = Satellite(orbit.position, orbit.velocity, central_body)

        return OrbitDataAccess(satellite)

# TODO: Write tests as I go.
# TODO: Add variable presets (Earth - ISS, Earth - Geostationary, Mars - Phobos etc).
# TODO: Write proper docstrings
if __name__ == "__main__":
    root = Tk()

    app: OrbitVisualiser = OrbitVisualiser(root)

    root.mainloop()