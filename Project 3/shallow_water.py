"""Optional task 3: The Shallow Water Equations"""
import sources
import numpy as np
import matplotlib.pyplot as plt

LENGTH = 5  # Length of the pool (in metres)
DX = 0.1  # (metres)
DT = 0.001  # (seconds)
T = np.arange(DT, 2 + DT, DT)  # NumPy array containing all discrete time steps
X = np.arange(0, LENGTH + DX, DX)  # NumPy array containing all discrete space steps
BETA = 1
G = 1


def h_u_derivatives(h, u, zeta, dx=DX, beta=BETA, g=G):
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
    def dy_dx(var):
        return (var[2:] - var[:-2]) / dx

    du_dx = dy_dx(u)
    dz_dx = dy_dx(zeta)

    dh_dt = dy_dx(h * u)
    du_dt = -du_dx * u - g * dz_dx - beta * u**3 / h**2
    return np.pad(dh_dt, 1, mode='edge'), np.pad(du_dt, 1)


if __name__ == '__main__':
    pass
