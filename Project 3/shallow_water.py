"""Optional task 3: The Shallow Water Equations"""
import sources
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

def gaussian(x, half_length):
    """Analytical single gaussian solution for advection in one dimension.

    Args:
        x: Array consisting of x.
        half_length: Half of the length L of the domain.

    Returns:
        A NumPy array of values of u along the x axis.
    """
    return np.exp(-(x - half_length)**2)

now = datetime.now()  # Get the current date and time
timestamp = now.strftime("%Y-%m-%d_%H.%M.%S")
LENGTH = 5  # Length of the pool (in metres)
DX = 0.1  # (metres)
DT = 0.001  # (seconds)
T = np.arange(DT, 2 + DT, DT)  # NumPy array containing all discrete time steps
X = np.arange(0, LENGTH + DX, DX)  # NumPy array containing all discrete space steps
BETA = 1
G = 1
U0 = np.zeros(len(X))
H0 = np.ones(len(X)) + gaussian(X, max(X)/2) #Due to considering B=1 we can say that dZETA/dx = dH/dx thus Zeta is not really neccessary anymore.
state0 = np.array([H0, U0])




def h_u_derivatives(state, dx=DX, beta=BETA, g=G):
    """Calculates the time derivatives of H and u.

    The derivative of u is padded at the boundaries with 0, since u=0 at the edges.
    The derivative of H is padded with the outer values, meaning the edge points move up and down
    in tandem with their neighbours.

    Args:
        h: Array of current depth at every spacial grid point.
        u: Array of current water velocity at every spatial grid point.
        zeta: Array of eviations from equilibrium surface height.
        dx: Spacial step size.
        beta: Constant value which dictates bottom friction.
        g: Gravitational accelleration constant.
    Returns:
        NumPy arrays of the derivative of H and u w.r.t. time at each x.
    """
    h, u = state
    zeta = state[0] - 1
    def dy_dx(var):
        return (var[2:] - var[:-2]) / dx

    du_dx = dy_dx(u)
    dz_dx = dy_dx(zeta)

    dh_dt = -dy_dx(h * u)
    du_dt = -du_dx * u - g * dz_dx - beta * u**3 / h**2
    return np.array([np.pad(dh_dt, 1, mode='edge'), np.pad(du_dt, 1)])

def numerical_data(method, state=state0, derivative=h_u_derivatives, dt=DT, t=T, dx=DX, beta=BETA, g=G):

    print("Working...\r", end='')
    data = np.empty((len(t), 2, len(state[0])), dtype=float)
    
    if method == sources.euler:
        for i in range(len(t)):
            data[i,:,:] = state
            state = method(state, t, derivative)
    else:
        print("this method has not been properly implemented yet.")
    return(data)

DATA_EULER = numerical_data(sources.euler)



if __name__ == '__main__':
    pass
