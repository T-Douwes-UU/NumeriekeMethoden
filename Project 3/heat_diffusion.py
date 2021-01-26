
"""Task 1: Heat diffusion"""
import sources
from sources import *
import numpy as np
from scipy.special import erfc
import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation

LENGTH = 5  # Length to plot (in metres), actual length is infinite
DX = 0.001 # (metres)
DT = 0.01 # (seconds)
T = np.arange(DT, 2 + DT, DT) # NumPy array containing all discrete time steps
X = np.arange(0, LENGTH + DX, DX) #NumPy array containing all discrete space steps
TEMP_0 = 0 # Initial temperature of the rod
TEMP_1 = 1 # The temperature at which the end of the rod is held at.
TEMP = np.full(len(X), TEMP_0, float)
TEMP[0] = TEMP_1
print(TEMP)
KAPPA = 1. #The thermal diffusion constant


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





def derivative(temp=TEMP, dt=DT, kappa=KAPPA):
    """
    Args:
        Current temperature at given time
    Return:
        A NumPy array of the double derivative of temperature to space 
    """
    new_var = (temp[2:] - 2 * temp[1:-1] + temp[:-2])/dt**2
    temp[1:-1] += new_var * KAPPA * dt
    temp[-1] = temp[-2]
    return temp
    

def produce_data(temp=TEMP, derivative=derivative, dt=DT, t=T):
    print("Working...")
    temp_list = np.zeros(,dtype=float)
    for i in T:
        new_temp = euler(temp, dt, derivative)
        temp = np.concatenate(temp, new_temp)
    print("Finished creating data array.")
    return temp

final_temp = produce_data()





def animate(length=LENGTH, temp_0=TEMP_0, temp_1=TEMP_1, kappa=KAPPA):
    """Returns a FuncAnimation object."""
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
    
    return sources.Player(fig, update, frames=len(T), interval=20)


if __name__ == '__main__':
    anim = animate()
    
    
    
