"""Task 2: Advection in 1D"""
import sources
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

LENGTH = 10  # Length to plot, actual length is infinite
WIDTH = 1  # Width of the Molenkamp solution 
DX = 0.01
DT = 0.001
X = np.arange(0, LENGTH + DX, DX)
T = np.arange(DT, 10+DT, DT)
CONST = 0.5  # traveling speed of variations in u.
now = datetime.now()  # Get the current date and time
timestamp = now.strftime("%Y-%m-%d_%H.%M.%S")

def gaussian(x_t, half_length):
    """Analytical single gaussian solution for advection in one dimension.

    Args:
        x_t: Array consisting of x - c * t (mod L).
        half_length: Half of the length L of the periodic domain.

    Returns:
        A NumPy array of values of u along the x axis.
    """
    return np.exp(-(x_t - half_length)**2)


def molenkamp(x_t, width):
    """Analytical Molenkamp solution for advection in one dimension.
    Args:
        x_t: Array consisting of x - c * t (mod L).
        width: Width of the triangle peak.
    Returns:
        A NumPy array of values of u along the x axis.
    """
    peak = np.minimum(x_t, width - x_t)
    return 2 * np.maximum(peak, 0) / width


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



def analytical_data(solution, t=T, x=X, length=LENGTH, width=WIDTH):
    if solution == gaussian:
        var = length/2
    elif solution == molenkamp:
        var = width
    data = np.array([solution(x-i, var) for i in t])
    return(data)


def numerical_data(method, state, derivative=u_derivative, dt=DT, t=T, dx=DX, const=CONST):
    """Constructs a 2D array containing numerically obtained data using a specified method.

    Args:
        state: Current state of u at every spacial grid point.
        derivative: spacial derivative of u.
        method: Chosen numerical method.
        dt: Time step size.
        t: Array containing all time steps.
        dx: Spatial step size.
        const: traveling speed of variations in u.
    Returns:
        A 2D NumPy array containing the state of u at every x step and every t step
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

#  Creating Data arrays of the numerical solutions.
DATA_GAUSS_ANLYTC = analytical_data(gaussian)
DATA_MOLEN_ANLYTC = analytical_data(molenkamp)

#  Creating data arrays using numerical methods on the Gaussian solution.
DATA_GAUSS_EULER = numerical_data(sources.euler, DATA_GAUSS_ANLYTC[0])
DATA_GAUSS_RK = numerical_data(sources.runge_kutta, DATA_GAUSS_ANLYTC[0])
DATA_GAUSS_LEAP = numerical_data(sources.leap_frog, DATA_GAUSS_ANLYTC[0])
DATA_GAUSS_ADAMS = numerical_data(sources.adams_bashforth, DATA_GAUSS_ANLYTC[0])
DATA_GAUSS_CRANK = numerical_data(crank_nicolson, DATA_GAUSS_ANLYTC[0])

#  Creating data arrays using numerical methods on the Molenkamp solution.
DATA_MOLEN_EULER = numerical_data(sources.euler, DATA_MOLEN_ANLYTC[0])
DATA_MOLEN_RK = numerical_data(sources.runge_kutta, DATA_MOLEN_ANLYTC[0])
DATA_MOLEN_LEAP = numerical_data(sources.leap_frog, DATA_MOLEN_ANLYTC[0])
DATA_MOLEN_ADAMS = numerical_data(sources.adams_bashforth, DATA_MOLEN_ANLYTC[0])
DATA_MOLEN_CRANK = numerical_data(crank_nicolson, DATA_MOLEN_ANLYTC[0])

# Calculate the <R^2> between the numerical methods and the Gaussian solution.
DEV_GAUSS_EULER = np.mean((DATA_GAUSS_EULER - DATA_GAUSS_ANLYTC)**2, axis=1)
DEV_GAUSS_RK = np.mean((DATA_GAUSS_RK - DATA_GAUSS_ANLYTC)**2, axis=1)
DEV_GAUSS_LEAP = np.mean((DATA_GAUSS_LEAP - DATA_GAUSS_ANLYTC)**2, axis=1)
DEV_GAUSS_ADAMS = np.mean((DATA_GAUSS_ADAMS - DATA_GAUSS_ANLYTC)**2, axis=1)
DEV_GAUSS_CRANK = np.mean((DATA_GAUSS_CRANK - DATA_GAUSS_ANLYTC)**2, axis=1)

# Calculate the <R^2> between the numerical methods and the Gaussian solution.
DEV_MOLEN_EULER = np.mean((DATA_MOLEN_EULER - DATA_MOLEN_ANLYTC)**2, axis=1)
DEV_MOLEN_RK = np.mean((DATA_MOLEN_RK - DATA_MOLEN_ANLYTC)**2, axis=1)
DEV_MOLEN_LEAP = np.mean((DATA_MOLEN_LEAP - DATA_MOLEN_ANLYTC)**2, axis=1)
DEV_MOLEN_ADAMS = np.mean((DATA_MOLEN_ADAMS - DATA_MOLEN_ANLYTC)**2, axis=1)
DEV_MOLEN_CRANK = np.mean((DATA_MOLEN_CRANK - DATA_MOLEN_ANLYTC)**2, axis=1)


def deviation(x=X, t=T):
    plt.figure()
    plt.yscale('log')
    plt.title(f"deviation between numerical and the Gaussian solution $<R^{2}>$")
    plt.xlabel("$t$ (s)")
    plt.ylabel("$<R^{2}>$")
    plt.plot(t, DEV_GAUSS_EULER, label='Euler method')
    plt.plot(t, DEV_GAUSS_RK, label='Runge-Kutta method')
    plt.plot(t, DEV_GAUSS_LEAP, label='Leap-frog method')
    plt.plot(t, DEV_GAUSS_ADAMS, label='Adams-Bashforth method')
    plt.plot(t, DEV_GAUSS_CRANK, label='Crank-Nicolson method')
    plt.legend()
    plt.savefig(f"deviation task2 Gaussian {timestamp}.png", dpi=300, bbox_inches='tight')
    plt.show()

    plt.figure()
    plt.yscale('log')
    plt.title(f"deviation between numerical and the Molenkamp solution $<R^{2}>$")
    plt.xlabel("$t$ (s)")
    plt.ylabel("$<R^{2}>$")
    plt.plot(t, DEV_MOLEN_EULER, label='Euler method', linestyle='--')
    plt.plot(t, DEV_MOLEN_RK, label='Runge-Kutta method')
    plt.plot(t, DEV_MOLEN_LEAP, label='Leap-frog method')
    plt.plot(t, DEV_MOLEN_ADAMS, label='Adams-Bashforth method')
    plt.plot(t, DEV_MOLEN_CRANK, label='Crank-Nicolson method')    
    plt.legend()
    plt.savefig(f"deviation task2 Molenkamp {timestamp}.png", dpi=300, bbox_inches='tight')
    plt.show()


def animate(length=LENGTH, width=WIDTH, x=X, t=T, const=CONST, beginframe=0):
    """Returns a FuncAnimation object."""
    #  x = np.arange(0, length, 0.01)
    #  t = np.arange(0.01, 10, 0.02)

    fig = plt.figure()
    ax = plt.axes(xlim=(0, length), ylim=(-0.2, 1.2))

    gauss_euler, = ax.plot([], [], label='Gaussian solution with Euler')
    gauss_rk, = ax.plot([], [], label='Gaussian solution with Runge-Kutta')
    gauss_leap, = ax.plot([], [], label='Gaussian solution with Leap-frog')
    gauss_adams, = ax.plot([], [], label='Gaussian solution with Adams-Bashforth')
    gauss_crank, = ax.plot([], [], c='violet', label='Gaussian solution with Crank-Nicolson')
    gauss, = ax.plot([], [], c='lime', linestyle=(0, (5, 5)), label='Analytical Gaussian solution')

    #molen_euler, = ax.plot([], [], label='Molenkamp solution with Euler')
    #molen_rk, = ax.plot([], [], label='Molenkamp solution with Runge-Kutta')
    #molen_leap, = ax.plot([], [], label='Molenkamp solution with Leap-frog')
    #molen_adams, = ax.plot([], [], label='Molenkamp solution with Adams-Bashforth')
    #molen_crank, = ax.plot([], [], c='violet', label='Molenkamp solution with Crank-Nicolson')
    #molen, = ax.plot([], [], c='lime', linestyle=(0, (5, 5)), label='Analytical Molenkamp solution')

    plt.legend(loc='upper left', fontsize='small')
    ax.set_title("Advection in one dimension")
    ax.set_xlabel("$x$ (m)")
    ax.set_ylabel("$u$")

    half_length = length / 2

    def update(i):
        x_t = (x - const * t[i]) % length  # Periodic domain

        gauss.set_data(x, gaussian(x_t, half_length))  # Update the plot
        gauss_euler.set_data(x, DATA_GAUSS_EULER[i])
        gauss_rk.set_data(x, DATA_GAUSS_RK[i])
        gauss_leap.set_data(x, DATA_GAUSS_LEAP[i])
        gauss_adams.set_data(x, DATA_GAUSS_ADAMS[i])
        gauss_crank.set_data(x, DATA_GAUSS_CRANK[i])

        #molen.set_data(x, molenkamp(x_t, width))
        #molen_euler.set_data(x, DATA_MOLEN_EULER[i])
        #molen_rk.set_data(x, DATA_MOLEN_RK[i])
        #molen_leap.set_data(x, DATA_MOLEN_LEAP[i])
        #molen_adams.set_data(x, DATA_MOLEN_ADAMS[i])
        #molen_crank.set_data(x, DATA_MOLEN_CRANK[i])
        
        #if you want to make a screenshot of a specific frame you can uncomment the following lines:
        #if i == 3000: #specify here what frame you would like to save.
        #    fig.savefig(f'plot task2 frame{i} {timestamp}.png', dpi=300, bbox_inches='tight')

        return gauss, gauss_euler, gauss_rk, gauss_leap, gauss_adams, gauss_crank, \
            molen, molen_euler, molen_rk, molen_leap, molen_adams, molen_crank

    #return sources.Player(fig, update, frames=len(t), interval=20) #Interactive animation for finding which frame you want.
    return FuncAnimation(fig, update, frames=range(beginframe,len(t)), interval=20, blit=True)


if __name__ == '__main__':
    anim = animate()

deviation()
