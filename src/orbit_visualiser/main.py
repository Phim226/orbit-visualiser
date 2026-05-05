import sys
from ttkbootstrap import Window
from orbit_visualiser.core import Orbit, Satellite, CentralBody
from orbit_visualiser.ui import UIController, UIBuilder, OrbitDataAccess, GeometryManager, initial_config

# TODO: Update propagation code and README description
class OrbitVisualiser():

    def __init__(self, root: Window, geo_manager: GeometryManager):
        oda: OrbitDataAccess = self._initialise_orbit_objects()

        builder = UIBuilder(root, oda, geo_manager)
        controller = UIController(builder, oda)
        builder.build(
            controller.reset_state,
            controller.validate_manual_input,
            controller.slider_changed,
            controller.format_display_value
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
    root = Window(title = "3D Orbit Visualiser", themename = "darkly", position = (0, 0))

    geo_manager = GeometryManager(sys.platform, root)

    app: OrbitVisualiser = OrbitVisualiser(root, geo_manager)

    root.mainloop()