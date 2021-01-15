"""Task 1: Heat diffusion"""
import sources
import numpy as np
from scipy.special import erfc
import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation


def analytical(x: np.ndarray, t: float, temp_0=0, temp_1=1, kappa=1):
    return temp_0 + (temp_1 - temp_0) * erfc(x / 2 / np.sqrt(kappa * t))


def animate():
    x = np.linspace(0, 5, 100)
    t = np.arange(0.01, 2, 0.01)

    fig = plt.figure()
    ax = plt.axes(xlim=(0, 5), ylim=(0, 1))
    line, = ax.plot([], [])
    ax.set_title("Heat diffusion in a rod")

    # animation function.  This is called sequentially
    def update(i):
        y = analytical(x, t[i])
        line.set_data(x, y)
        return line,

    return sources.Player(fig, update, frames=len(t), interval=20)


if __name__ == '__main__':
    anim = animate()
