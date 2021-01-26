
"""Task 1: Heat diffusion"""
import sources
from sources import *
import numpy as np
from scipy.special import erfc
import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation

LENGTH = 5  # Length to plot (in metres), actual length is infinite
DX = 0.1  # (metres)
DT = 0.001  # (seconds)
T = np.arange(DT, 2 + DT, DT)  # NumPy array containing all discrete time steps
X = np.arange(0, LENGTH + DX, DX)  # NumPy array containing all discrete space steps
TEMP_0 = 0  # Initial temperature of the rod
TEMP_1 = 1  # The temperature at which the end of the rod is held at.
TEMP = np.full(len(X), TEMP_0, float)
TEMP[0] = TEMP_1
print(TEMP)
KAPPA = 1.  # The thermal diffusion constant


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


def derivative(temp, dt=DT, dx=DX, kappa=KAPPA):
    """
    Args:
        Current temperature at given time
    Return:
        A NumPy array of the double derivative of temperature to space 
    """
    d_temp = kappa * (temp[2:] - 2 * temp[1:-1] + temp[:-2]) / dx**2
    return np.pad(d_temp, 1, constant_values=(0, d_temp[-1]))


def produce_data(temp=TEMP, derivative=derivative, dt=DT, t=T):
    print("Working...")
    temp_list = np.zeros((len(t), len(temp)), dtype=float)
    for i in range(len(t)):
        temp_list[i] = temp
        temp = euler(temp, dt, derivative)

    print("Finished creating data array.")
    return temp_list


final_temp = produce_data()


def animate(x=X, t=T, length=LENGTH, temp_0=TEMP_0, temp_1=TEMP_1, kappa=KAPPA):
    """Returns a FuncAnimation object."""
    fig = plt.figure()
    ax = plt.axes(xlim=(0, length), ylim=(temp_0, temp_1))
    anlytc, = ax.plot([], [])
    eul, = ax.plot([], [])
    ax.set_title("Heat diffusion in a half-infinite rod")
    ax.set_xlabel("$x$ (m)")
    ax.set_ylabel("$T$ (${^\circ}$C)")

    def update(i):
        y = analytical(x, t[i], temp_0, temp_1, kappa)  # Calculate the analytical temperature
        anlytc.set_data(x, y)  # Update the plot
        eul.set_data(x, final_temp[i])
        return anlytc,
    
    return sources.Player(fig, update, frames=len(T), interval=20)


if __name__ == '__main__':
    anim = animate()
