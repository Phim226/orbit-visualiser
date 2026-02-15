import sys
from tkinter import Tk
from orbit_visualiser.core import Orbit, Satellite
from orbit_visualiser.core.neworbit import NewOrbit, CentralBody
from orbit_visualiser.core.satellite import NewSatellite
from orbit_visualiser.ui import OrbitFigure, OrbitConfigBuilder, OrbitConfigController
from orbit_visualiser.ui.common.presets import initial_config

class OrbitVisualiser():

    figure_frame_placement = ("left", "nw")
    config_frame_placement = ("right", "ne")

    def __init__(self, root: Tk):
        root.title("2D Orbit Visualiser")

        if sys.platform.startswith("win"):
            root.state("zoomed")
        else:
            root.state("normal")

        orbit: NewOrbit = NewOrbit.from_orbital_elements(
            initial_config.eccentricity,
            initial_config.radius_of_periapsis,
            initial_config.gravitational_parameter,
            initial_config.true_anomaly,
        )

        central_body: CentralBody = CentralBody(
            initial_config.gravitational_parameter,
            initial_config.radius
        )
        satellite: NewSatellite = NewSatellite(orbit.position, orbit.velocity, central_body)

        orbit_figure: OrbitFigure = OrbitFigure(
            root, self.figure_frame_placement, satellite
        )
        orbit_figure.build()

        orbit_builder: OrbitConfigBuilder = OrbitConfigBuilder(
            root, self.config_frame_placement, orbit, central_body, satellite
        )
        orbit_controller: OrbitConfigController = OrbitConfigController(
            orbit_figure, orbit_builder, orbit, satellite, central_body
        )
        orbit_builder.build(
            orbit_controller.reset_state,
            orbit_controller.validate_manual_input,
            orbit_controller.slider_changed,
            orbit_controller.format_display_value
        )

# TODO: Write tests as I go.
# TODO: Add variable presets (Earth - ISS, Earth - Geostationary, Mars - Phobos etc).
# TODO: Write proper docstrings
if __name__ == "__main__":
    root = Tk()

    app: OrbitVisualiser = OrbitVisualiser(root)

    root.mainloop()