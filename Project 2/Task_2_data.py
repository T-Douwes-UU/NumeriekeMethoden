"""
This file produces the data for Task 2. Running it prompts the user to input initial values for the system,
after which the necessary data is produced and output to a file in the parent directory.

This file must be in the same folder as sources.py to be run successfully.
"""
from sources import *

# ==================================================
# #### Defining constants #####
# ==================================================
G = 1
omega = 1  # angular velocity
m_1 = set_value("Please enter the mass of the first primary in kg: ")  # = 1.9885E30 (Sun)
m_2 = set_value("Please enter the mass of the second primary in kg: ")  # = 1.89819E27 (Jupiter)
m1 = m_1 / (m_1 + m_2)
m2 = m_2 / (m_1 + m_2)
END_TIME = 1000


# =========================================================
# #### Defining derivatives of the object's trajectory#####
# =========================================================
def derivatives_threebody(variables, *_, g=G):
    """Calculates the derivatives with respect to time of the variables in the restricted three body system.

    We make use of the EOM for a corotating frame of reference in which the two primary masses are stationary.
    Note that these EOM have no explicit time depencency; the second positional argument of this function is
    in place for compatibility with runge_kutta(), but we don't need any actual time input.

    Args:
        variables: NumPy array containing the relative variables x, y, vx, and vy
        g (optional): Gravitational constant G in Newton's law of gravitation.

    Returns:
        A NumPy array containing the derivatives dx/dt, dy/dt, dvx/dt, and dvy/dt.
    """
    x, y, vx, vy = variables

    # for use in dvx and dvy
    x1, x2, y1, y2 = positions(1, 0, m1, m2)
    r1, r2 = abs(x1), abs(x2)
    R1 = np.sqrt((x - x1)**2 + (y - y1)**2)
    R2 = np.sqrt((x - x2)**2 + (y - y2)**2)
    const1 = -g * m1 / R1**3   
    const2 = -g * m2 / R2**3

    # We calculate the actual derivatives
    dx = vx
    dy = vy
    dvx = const1 * (x - r1) + const2 * (x + r2) + omega**2 * x + 2*omega*vy
    dvy = y * (const1 + const2) + omega**2 * y - 2*omega*vx

    return np.array([dx, dy, dvx, dvy])


# ==================================================
# #### Producing data file and running simulation#####
# ==================================================
def main():
    x = set_value("Please enter the initial x position of the asteroid: ")
    y = set_value("Please enter the initial y position of the asteroid: ")
    vx = set_value("Please enter the initial x velocity of the asteroid: ")
    vy = set_value("Please enter the initial y velocity of the asteroid: ")
    dt = set_value("Please enter the time step dt for integration in years: ")
    variables = np.array([x, y, vx, vy])
    produce_data(2, variables, derivatives_threebody, END_TIME, dt, m1, m2)


main()
