
"""Task 1: Heat diffusion"""
import sources
from sources import *
import numpy as np
from scipy.special import erfc
import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation

LENGTH = 5  # Length to plot (in metres), actual length is infinite
DX = 0.1 # (metres)
DT = 0.001 # (seconds)
T = np.arange(DT, 2 + DT, DT) # NumPy array containing all discrete time steps
X = np.arange(0, LENGTH + DX, DX) #NumPy array containing all discrete space steps
TEMP_0 = 0 # Initial temperature of the rod
TEMP_1 = 1 # The temperature at which the end of the rod is held at.
TEMP = np.full(len(X), TEMP_0, float)
TEMP[0] = TEMP_1
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


def derivative(temp, dt=DT, dx=DX, kappa=KAPPA):
    """
    Args:
        temp: Current temperature at every spacial grid point.
        dt: Timestep size
        dx: Spacial step size
        kappa: The thermal diffusion constant.
    Return:
        A NumPy array of the double derivative of temperature to space 
    """
    d_temp = kappa * (temp[2:] - 2 * temp[1:-1] + temp[:-2]) / dx**2
    return np.pad(d_temp, 1, constant_values=(0, d_temp[-1]))
    

def numerical_data(temp=TEMP, derivative=derivative, method=euler, dt=DT, t=T):
    """
    Args:
        temp: Current temperature at every spacial grid point.
        derivative: Double spacial derivative of the temperature.
        method: Chosen numerical method.
        dt: Timestep size.
        t: Array containing all time steps.
    Return:
        A NumPy array containing the temperatures at every x step and every t step
        found using the chosen method method.
    """
    print("Working...")
    temp_list = np.empty((len(t), len(temp)), dtype=float)


    if method in (leap_frog, adams_bashforth):
        temp_list[0:2,:] = kutta_begin(temp, dt, derivative)
        temp = temp_list[1]
        #t = t[1:]
        
        for i in range(len(t[:-2])):
            temp = method(temp_list[i+1], temp_list[i], dt, derivative)
            temp_list[i+2] = temp
        
    else:
        for i in range(len(t)):    
            temp_list[i] = temp
            temp = method(temp, dt, derivative)
        
    print("Finished creating data array.")
    return temp_list

TEMP_euler = numerical_data(method=euler)
TEMP_RK = numerical_data(method=runge_kutta)
TEMP_LEAP = numerical_data(method=leap_frog)

def animate(x=X, t=T, length=LENGTH, temp_0=TEMP_0, temp_1=TEMP_1, kappa=KAPPA):
    """Returns a FuncAnimation object."""
    fig = plt.figure()
    ax = plt.axes(xlim=(0, length), ylim=(temp_0, temp_1))
    eul, = ax.plot([], [], label='euler method')
    RK, = ax.plot([], [], label='Runge-Kutta method')
    LEAP, = ax.plot([], [], label='Leap-Frog method')
    anlytc, = ax.plot([], [], label='anlytical result', linestyle='dashed')
    plt.legend()
    ax.set_title("Heat diffusion in a half-infinite rod")
    ax.set_xlabel("$x$ (m)")
    ax.set_ylabel("$T$ (${^\circ}$C)")

    def update(i):
        y = analytical(x, t[i], temp_0, temp_1, kappa)  # Calculate the analytical temperature
        anlytc.set_data(x, y)  # Update the plot
        eul.set_data(x, TEMP_euler[i])
        RK.set_data(x, TEMP_RK[i])
        LEAP.set_data(x, TEMP_LEAP[i])
        return anlytc,
    
    return sources.Player(fig, update, frames=len(T), interval=20)


if __name__ == '__main__':
    anim = animate()

