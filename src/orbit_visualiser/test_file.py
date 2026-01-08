import tkinter
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from math import pi

rp = 2 # radius of periapsis
ra = 8 # radius of apoapsis
e = 0.6 # eccentricity
a=5  # semimajor axis
b=4  # semiminor axis
p = (b**2)/(3*a - ra - rp) # orbital parameter
t = np.linspace(0, 2*pi, 100)

root = tkinter.Tk()
root.wm_title("Embedded in Tk")

fig = Figure(figsize=(5, 4), dpi=100)
fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
ax = fig.add_subplot()
ax.spines['left'].set_position(('data', 0))
ax.spines['bottom'].set_position(('data', 0))

ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')

ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
line, = ax.plot(p*(np.cos(t)/(1+e*np.cos(t))) , p*(np.sin(t)/(1+e*np.cos(t))))

canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
canvas.draw()

# pack_toolbar=False will make it easier to use a layout manager later on.
toolbar = NavigationToolbar2Tk(canvas, root, pack_toolbar=False)
toolbar.update()


def update_eccentricity(new_val):
    # retrieve frequency
    f = float(new_val)

    # update data
    x = p*(np.cos(t)/(1+f*np.cos(t)))
    y = p*(np.sin(t)/(1+f*np.cos(t)))
    line.set_data(x, y)

    # required to update canvas and attached toolbar!
    canvas.draw()


slider_update = tkinter.Scale(root, from_ = 0, to = 1, resolution = 0.1, orient = tkinter.HORIZONTAL,
                              command = update_eccentricity, label = "Eccentricity")
slider_update.set(0.6)

# Packing order is important. Widgets are processed sequentially and if there
# is no space left, because the window is too small, they are not displayed.
# The canvas is rather flexible in its size, so we pack it last which makes
# sure the UI controls are displayed as long as possible.
slider_update.pack(side=tkinter.RIGHT)
toolbar.pack(side=tkinter.BOTTOM, fill=tkinter.X)
canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=True)

tkinter.mainloop()