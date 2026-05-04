import sys
from tkinter import Tk
from orbit_visualiser.core import Orbit, Satellite, CentralBody
from orbit_visualiser.ui import UIController, UIBuilder, OrbitDataAccess, initial_config

# TODO: Fix broken tests

class OrbitVisualiser():

    INPUT_GEOMETRY = ("left", "nw")
    FIGURE_GEOMETRY = ("left", "nw")
    PROPS_GEOMETRY = ("left", "nw")

    def __init__(self, root: Tk):
        root.title("2D Orbit Visualiser")

        if sys.platform.startswith("win"):
            root.state("zoomed")
        else:
            root.state("normal")

        oda: OrbitDataAccess = self._initialise_orbit_objects()

        builder = UIBuilder(root, oda)
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
    root = Tk()

    app: OrbitVisualiser = OrbitVisualiser(root)

    root.mainloop()