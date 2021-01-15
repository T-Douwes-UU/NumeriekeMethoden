"""Task 1: Heat diffusion"""
import sources
import numpy as np
from scipy.special import erfc
import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation

LENGTH = 5  # Length to plot, actual length is infinite
T_0 = 0
T_1 = 1
KAPPA = 1


def analytical(x: np.ndarray, t, temp_0, temp_1, kappa):
    """Analytical solution for heat diffusion in a semi-infinite rod.

    Args:
        x: An array of distances x from the endpoint of the rod.
        t: The current time.
        temp_0: The initial temperature of the rod.
        temp_1: The temperature the end of the rod is held at.
        kappa: The thermal diffusion constant.

    Returns:
        A NumPy array of current temperatures along the x axis.
    """
    return temp_0 + (temp_1 - temp_0) * erfc(x / 2 / np.sqrt(kappa * t))


def animate(length=LENGTH, temp_0=T_0, temp_1=T_1, kappa=KAPPA):
    """Returns a FuncAnimation object."""
    x = np.linspace(0, length, 100)
    t = np.arange(0.01, 2, 0.01)

    fig = plt.figure()
    ax = plt.axes(xlim=(0, length), ylim=(temp_0, temp_1))
    line, = ax.plot([], [])
    ax.set_title("Heat diffusion in a half-infinite rod")
    ax.set_xlabel("$x$ (m)")
    ax.set_ylabel("$T$ (${^\circ}$C)")

    def update(i):
        y = analytical(x, t[i], temp_0, temp_1, kappa)  # Calculate the analytical temperature
        line.set_data(x, y)  # Update the plot
        return line,

    return sources.Player(fig, update, frames=len(t), interval=20)


if __name__ == '__main__':
    anim = animate()
