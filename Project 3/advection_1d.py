"""Task 2: Advection in 1D"""
import sources
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

LENGTH = 10  # Length to plot, actual length is infinite
WIDTH = 1  # Width of the Molenkamp solution


def gaussian(x_t, half_length):
    """Analytical single gaussian solution for advection in one dimension.

    Args:
        x_t: Array consisting of x - t (mod L).
        half_length: Half of the length L of the periodic domain.

    Returns:
        A NumPy array of values of u along the x axis.
    """
    return np.exp(-(x_t - half_length)**2)


def molenkamp(x_t, width):
    """Analytical Molenkamp solution for advection in one dimension.

    Args:
        x_t: Array consisting of x - t (mod L).
        width: Width of the triangle peak.

    Returns:
        A NumPy array of values of u along the x axis.
    """
    return 2 * np.clip(x_t, 0, width - x_t) / width


def animate(length=LENGTH, width=WIDTH):
    """Returns a FuncAnimation object."""
    x = np.arange(0, length, 0.01)
    t = np.arange(0.01, 10, 0.02)

    fig = plt.figure()
    ax = plt.axes(xlim=(0, length), ylim=(0, 1.2))
    gauss, = ax.plot([], [])
    molen, = ax.plot([], [])
    ax.set_title("Advection in one dimension")
    ax.set_xlabel("$x$ (m)")
    ax.set_ylabel("$u$")

    half_length = length / 2

    def update(i):
        x_t = (x - t[i]) % length  # Periodic domain
        gauss.set_data(x, gaussian(x_t, half_length))  # Update the plot
        molen.set_data(x, molenkamp(x_t, width))
        return gauss, molen

    return FuncAnimation(fig, update, frames=len(t), interval=20, blit=True)


if __name__ == '__main__':
    anim = animate()
