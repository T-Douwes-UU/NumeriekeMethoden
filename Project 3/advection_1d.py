"""Task 2: Advection in 1D"""
import sources
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

LENGTH = 10  # Length to plot, actual length is infinite
WIDTH = 1  # Width of the Molenkamp solution 
DX = 0.01
DT = 0.001
X = np.arange(0, LENGTH + DX, DX)
T = np.arange(DT, 10+DT, DT)
CONST = 1  # traveling speed of variations in u.


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


def u_derivative(u, dx=DX, c=CONST):
    """Time derivative of u.
    
    Args:
        u: Array containing the current displacement u at every point x.
        dx: Spatial step size.
        c: Propagation speed.
    Returns:
        A NumPy array containing values for du/dt for each point x.
    """
    return -c * (np.roll(u, -1) - np.roll(u, 1)) / (2 * dx)


def crank_nicolson(state, dt=DT, const=CONST, dx=DX):
    """Constructs the matrix C used in the Crank-Nicolson integration method.

    Args:
        state: Array of values from which the length is extracted.
        dt: Time step.
        const: Constant coefficient in the PDE.
        dx: Spatial step.
    Returns:
        A 2D NumPy array representing the matrix C.
    """
    c = const * dt / (4 * dx)
    n = len(state)
    i = np.identity(n)

    #  Construct A and B from diagonal matrices
    arr = c * np.eye(n, k=1)
    array = arr - arr.T
    array[-1, 0] = c
    array[0, -1] = -c
    
    A = i + array
    B = i - array
    C = np.linalg.inv(A) @ B
    return C


STATE_GAUSS = gaussian(X, LENGTH/2)
STATE_MOLEN = molenkamp(X, WIDTH)


def numerical_data(method, state, derivative=u_derivative, dt=DT, t=T, dx=DX, const=CONST):
    """Constructs a 2D array containing numerically obtained data using a specified method.

    Args:
        state: Current temperature at every spacial grid point.
        derivative: Double spacial derivative of the temperature.
        method: Chosen numerical method.
        dt: Time step size.
        t: Array containing all time steps.
        dx: Spatial step size.
        const: traveling speed of variations in u.
    Returns:
        A 2D NumPy array containing the temperatures at every x step and every t step
        found using the chosen method method.
    """
    print("Working...\r", end='')
    data = np.empty((len(t), len(state)), dtype=float)

    if method in (sources.leap_frog, sources.adams_bashforth):
        data[0] = state
        data[1] = sources.euler(state, dt, derivative)

        for i in range(len(t) - 2):
            data[i+2] = method(data[i+1], data[i], dt, derivative)

    elif method == crank_nicolson:
        C = crank_nicolson(state, dt, const, dx)
        for i in range(len(t)):
            data[i] = state
            state = C.dot(state)

    else:
        for i in range(len(t)):
            data[i] = state
            state = method(state, dt, derivative)

    print(f"Finished creating data array using {method.__name__}.")
    return data


DATA_GAUSS_EULER = numerical_data(sources.euler, STATE_GAUSS)
DATA_GAUSS_RK = numerical_data(sources.runge_kutta, STATE_GAUSS)
DATA_GAUSS_LEAP = numerical_data(sources.leap_frog, STATE_GAUSS)
DATA_GAUSS_ADAMS = numerical_data(sources.adams_bashforth, STATE_GAUSS)
DATA_GAUSS_CRANK = numerical_data(crank_nicolson, STATE_GAUSS)

DATA_MOLEN_EULER = numerical_data(sources.euler, STATE_MOLEN)
DATA_MOLEN_RK = numerical_data(sources.runge_kutta, STATE_MOLEN)
DATA_MOLEN_LEAP = numerical_data(sources.leap_frog, STATE_MOLEN)
DATA_MOLEN_ADAMS = numerical_data(sources.adams_bashforth, STATE_MOLEN)
DATA_MOLEN_CRANK = numerical_data(crank_nicolson, STATE_MOLEN)


def animate(length=LENGTH, width=WIDTH, x=X, t=T, const=CONST):
    """Returns a FuncAnimation object."""
    #  x = np.arange(0, length, 0.01)
    #  t = np.arange(0.01, 10, 0.02)

    fig = plt.figure()
    ax = plt.axes(xlim=(0, length), ylim=(0, 1.2))
    gauss, = ax.plot([], [])
    gauss_euler, = ax.plot([], [])
    molen, = ax.plot([], [])
    molen_euler, = ax.plot([], [])
    ax.set_title("Advection in one dimension")
    ax.set_xlabel("$x$ (m)")
    ax.set_ylabel("$u$")

    half_length = length / 2

    def update(i):
        x_t = (x - const * t[i]) % length  # Periodic domain
        gauss.set_data(x, gaussian(x_t, half_length))  # Update the plot
        gauss_euler.set_data(x, DATA_GAUSS_EULER[i])
        molen.set_data(x, molenkamp(x_t, width))
        molen_euler.set_data(x, DATA_MOLEN_EULER[i])
        return gauss, gauss_euler, molen, molen_euler

    return FuncAnimation(fig, update, frames=len(t), interval=20, blit=True)


if __name__ == '__main__':
    anim = animate()
