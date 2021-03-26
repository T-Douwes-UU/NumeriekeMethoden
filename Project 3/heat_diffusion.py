"""Task 1: Heat diffusion"""
import sources
from datetime import datetime
import numpy as np
from scipy.special import erfc
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

LENGTH = 5  # Length to plot (in metres), actual length is infinite
DX = 0.1  # (metres)
DT = 0.001  # (seconds)
T = np.arange(DT, 2 + DT, DT)  # NumPy array containing all discrete time steps
X = np.arange(0, LENGTH + DX, DX)  # NumPy array containing all discrete space steps
TEMP_0 = 0  # Initial temperature of the rod
TEMP_1 = 1  # The temperature at which the end of the rod is held at.
TEMP = np.full(len(X), TEMP_0, float)
TEMP[0] = TEMP_1
KAPPA = sources.set_value('Please enter the value for Kappa: ')  # The thermal diffusion constant
now = datetime.now()  # Get the current date and time
timestamp = now.strftime("%Y-%m-%d_%H.%M.%S")
stability = KAPPA * DT / DX**2


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
        A NumPy array of the derivative of temperature w.r.t. time at each x.
    """
    d_temp = kappa * (temp[2:] - 2 * temp[1:-1] + temp[:-2]) / dx**2
    return np.pad(d_temp, 1, constant_values=(0, d_temp[-1]))


def crank_nicolson(state, dt, const, dx):
    """Constructs the matrix C used in the Crank-Nicolson integration method.

    Args:
        state: Array of values from which the length is extracted.
        dt: Time step.
        const: Constant coefficient in the PDE.
        dx: Spatial step.
    Returns:
        A 2D NumPy array representing the matrix C.
    """
    c = const * dt / (2 * dx**2)
    n = len(state)
    i = np.identity(n)

    # Construct A and B from diagonal matrices
    arr1 = -2 * c * i
    arr2 = c * np.eye(n, k=1)
    array = arr1 + arr2 + arr2.T
    array[0] = 0
    array[-1] = array[-2]

    A = i - array
    B = i + array
    C = np.linalg.inv(A) @ (B)
    return C


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
        found using the chosen method method.
    """
    print("Working...\r", end='')
    data = np.empty((len(t), len(temp)), dtype=float)

    if method in (sources.leap_frog, sources.adams_bashforth):
        data[0] = temp
        data[1] = sources.euler(temp, dt, derivative)

        for i in range(len(t) - 2):
            data[i+2] = method(data[i+1], data[i], dt, derivative)

    elif method == crank_nicolson:
        C = method(temp, dt, kappa, dx)
        for i in range(len(t)):
            data[i] = temp
            temp = C.dot(temp)

    else:
        for i in range(len(t)):
            data[i] = temp
            temp = method(temp, dt, derivative)

    print(f"Finished creating data array using {method}.")
    return data

DATA_anlytc = np.array([analytical(X, i, TEMP_0, TEMP_1, KAPPA) for i in T]) #  Calculate the analytical temperature in all time steps.

# Calculate the temperatures using various numerical methods.
DATA_euler = numerical_data(sources.euler)
DATA_RK = numerical_data(sources.runge_kutta)
DATA_LEAP = numerical_data(sources.leap_frog)
DATA_ADAMS = numerical_data(sources.adams_bashforth)
DATA_CRANK = numerical_data(crank_nicolson)

# Calculate the <R^2> between the numerical methods and analytical solution.
DEV_euler = np.mean((DATA_euler - DATA_anlytc)**2, axis=1)
DEV_RK = np.mean((DATA_RK - DATA_anlytc)**2, axis=1)
DEV_LEAP = np.mean((DATA_LEAP - DATA_anlytc)**2, axis=1)
DEV_ADAMS = np.mean((DATA_ADAMS - DATA_anlytc)**2, axis=1)
DEV_CRANK = np.mean((DATA_CRANK - DATA_anlytc)**2, axis=1)


def deviation(x=X, t=T):
    plt.figure()
    plt.yscale('log')
    plt.title(rf"deviation between numerical and analytical solution $<R^{{2}}>$ with $(\frac{{\kappa \Delta t}}{{\Delta x^2}})$ = {stability:.3f}", wrap=True)
    plt.xlabel("$t$ (s)")
    plt.ylabel("$<R^{2}>$")
    plt.plot(t, DEV_euler, label='Euler method')
    plt.plot(t, DEV_RK, label='Runge-Kutta method')
    plt.plot(t, DEV_LEAP, label='Leap-frog method')
    plt.plot(t, DEV_ADAMS, label='Adams-Bashforth method')
    plt.plot(t, DEV_CRANK, label='Crank-Nicolson method')
    plt.legend(loc='lower right', fontsize='small')
    plt.savefig(f"Deviation plot task 1 {timestamp}.png", dpi=300, bbox_inches='tight')


def animate(x=X, t=T, length=LENGTH, temp_0=TEMP_0, temp_1=TEMP_1, kappa=KAPPA, beginframe=0):
    """Returns a FuncAnimation object."""
    fig = plt.figure()
    ax = plt.axes(xlim=(0, length), ylim=(temp_0, temp_1))
    eul, = ax.plot([], [], label='euler method')
    rk, = ax.plot([], [], label='Runge-Kutta method')
    leap, = ax.plot([], [], label='Leap-Frog method')
    adams, = ax.plot([], [], label='Adams-Bashforth method')
    crank, = ax.plot([], [], c='violet', label='Crank-Nicolson method')
    anlytc, = ax.plot([], [], c='cyan', linestyle=(0, (5, 2)), label='Analytical result')
    plt.legend()
    ax.set_title(rf"Heat diffusion in a half-infinite rod with $(\frac{{\kappa \Delta t}}{{\Delta x^2}})$ = {stability:.3f}")
    ax.set_xlabel("$x$ (m)")
    ax.set_ylabel("$T$ (${^\circ}$C)")

    def update(i):
        anlytc.set_data(x, DATA_anlytc[i])  # Update the plot
        eul.set_data(x, DATA_euler[i])
        rk.set_data(x, DATA_RK[i])
        leap.set_data(x, DATA_LEAP[i])
        adams.set_data(x, DATA_ADAMS[i])
        crank.set_data(x, DATA_CRANK[i])
        #if you want to make a screenshot of a specific frame you can uncomment the following lines:
        #if i == 60: #specify here what frame you would like to save.
        #    fig.savefig(f'plot frame{i} of task1 {timestamp}.png', dpi=300, bbox_inches='tight')

    return sources.Player(fig, update, frames=len(T), interval=20) #Interactive animation.
    #If you just want the plain animation, use this instead of the sources.player function:
    #return FuncAnimation(fig, update, frames=range(beginframe,len(T)), interval=20)


if __name__ == '__main__':
    anim = animate()

deviation()
