from tkinter import Tk
from orbit_visualiser.core import Orbit
from orbit_visualiser.ui import OrbitFigure
from orbit_visualiser.ui import OrbitConfigurer

class OrbitVisualiser():

    figure_frame_placement = ("left", "nw")
    config_frame_placement = ("right", "ne")

    def __init__(self, root: Tk):
        root.title("2D Orbit Visualiser")

        orbit: Orbit = Orbit()

        orbit_figure: OrbitFigure = OrbitFigure(root, self.figure_frame_placement, orbit)
        orbit_figure.build()

        orbit_config: OrbitConfigurer = OrbitConfigurer(root, self.config_frame_placement, orbit_figure, orbit)
        orbit_config.build()

# TODO: Write tests as I go.
if __name__ == "__main__":
    root = Tk()

    app: OrbitVisualiser = OrbitVisualiser(root)

    root.mainloop()