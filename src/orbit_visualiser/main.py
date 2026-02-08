from tkinter import Tk
from orbit_visualiser.core import Orbit, Satellite, CentralBody
from orbit_visualiser.ui import OrbitFigure, OrbitConfigBuilder, OrbitConfigController

class OrbitVisualiser():

    figure_frame_placement = ("left", "nw")
    config_frame_placement = ("right", "ne")

    def __init__(self, root: Tk):
        root.title("2D Orbit Visualiser")
        root.state("zoomed")

        orbit: Orbit = Orbit()
        central_body: CentralBody = CentralBody()
        satellite: Satellite = Satellite(orbit, central_body)

        orbit_figure: OrbitFigure = OrbitFigure(
            root, self.figure_frame_placement, orbit, central_body, satellite
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
            orbit_controller.update_value,
            orbit_controller.format_display_value
        )

# TODO: Write tests as I go.
# TODO: Add variable presets (Earth - ISS, Earth - Geostationary, Mars - Phobos etc).
if __name__ == "__main__":
    root = Tk()

    app: OrbitVisualiser = OrbitVisualiser(root)

    root.mainloop()