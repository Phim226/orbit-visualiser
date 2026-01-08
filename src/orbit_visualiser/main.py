from tkinter import Tk
from orbit_visualiser.core import Orbit
from orbit_visualiser.ui import OrbitFigure

class OrbitVisualiser():

    figure_frame_placement = ("left", "nw")
    config_frame_placement = ("right", "ne")

    def __init__(self, root: Tk):
        root.title("2D Orbit Visualiser")

        orbit: Orbit = Orbit()
        orbit_figure: OrbitFigure = OrbitFigure(root, self.figure_frame_placement, orbit)
        orbit_figure.build()


if __name__ == "__main__":
    root = Tk()

    app: OrbitVisualiser = OrbitVisualiser(root)

    root.mainloop()