"""Task 1: Heat diffusion"""
import sources
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


def temp_derivative(temp, dx=DX, kappa=KAPPA):
    """Calculates the time derivative of temperature, using the second spatial derivative.
    
    Args:
        temp: Current temperature at every spacial grid point.
        dx: Spacial step size.
        kappa: The thermal diffusion constant.
    Returns:
        A NumPy array of the derivative of temperature w.r.t. time at each x 
    """
    d_temp = kappa * (temp[2:] - 2 * temp[1:-1] + temp[:-2]) / dx**2
    return np.pad(d_temp, 1, constant_values=(0, d_temp[-1]))
    

def numerical_data(method, temp=TEMP, derivative=temp_derivative, dt=DT, t=T, dx=DX, kappa=KAPPA):
    """Constructs a 2D array containing numerically obtained data using a specified method.
    
    Args:
        temp: Current temperature at every spacial grid point.
        derivative: Double spacial derivative of the temperature.
        method: Chosen numerical method.
        dt: Time step size.
        t: Array containing all time steps.
        dx: Spatial step size.
        kappa: Thermal diffusion constant.
    Returns:
        A 2D NumPy array containing the temperatures at every x step and every t step
    """
    print("Working...")
    data = np.empty((len(t), len(temp)), dtype=float)

    if method in (sources.leap_frog, sources.adams_bashforth):
        data[0] = temp
        data[1] = sources.euler(temp, dt, derivative)

        for i in range(len(t)-2):
            data[i+2] = method(data[i+1], data[i], dt, derivative)
    
    elif method == sources.crank_nicolson:
        C = sources.crank_nicolson(temp, dt, derivative, kappa, dx)
        for i in range(len(t)):
            data[i] = temp
            temp = C.dot(temp)

    else:
        for i in range(len(t)):
            data[i] = temp
            temp = method(temp, dt, derivative)

    print(f"Finished creating data array using {method}.")
    return data


TEMP_euler = numerical_data(method=sources.euler)
#TEMP_RK = numerical_data(method=sources.runge_kutta)
#TEMP_LEAP = numerical_data(method=sources.leap_frog)
#TEMP_ADAMS = numerical_data(method=sources.adams_bashforth)
TEMP_CRANK = numerical_data(method=sources.crank_nicolson)


def animate(x=X, t=T, length=LENGTH, temp_0=TEMP_0, temp_1=TEMP_1, kappa=KAPPA):
    """Returns a FuncAnimation object."""
    fig = plt.figure()
    ax = plt.axes(xlim=(0, length), ylim=(temp_0, temp_1))
    eul, = ax.plot([], [], label='euler method')
    #  rk, = ax.plot([], [], label='Runge-Kutta method')
    #  leap, = ax.plot([], [], label='Leap-Frog method')
    #  adams, = ax.plot([],[], label='Adams-Bashforth method')
    crank, = ax.plot([],[], label='Crank-Nicolson method')
    anlytc, = ax.plot([], [], label='analytical result', linestyle='dashed')
    plt.legend()
    ax.set_title("Heat diffusion in a half-infinite rod")
    ax.set_xlabel("$x$ (m)")
    ax.set_ylabel("$T$ (${^\circ}$C)")

    def update(i):
        y = analytical(x, t[i], temp_0, temp_1, kappa)  # Calculate the analytical temperature
        anlytc.set_data(x, y)  # Update the plot
        eul.set_data(x, TEMP_euler[i])
        #  rk.set_data(x, TEMP_RK[i])
        #  leap.set_data(x, TEMP_LEAP[i])
        #  adams.set_data(x, TEMP_ADAMS[i])
        crank.set_data(x, TEMP_CRANK[i])
        return anlytc,
    
    return sources.Player(fig, update, frames=len(T), interval=20)


if __name__ == '__main__':
    anim = animate()
