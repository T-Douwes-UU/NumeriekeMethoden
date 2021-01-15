"""
This file produces the data for Task 1. Running it prompts the user to input initial values for the system,
after which the necessary data is produced and output to a file in the parent directory.

This file must be in the same folder as sources.py to be run successfully.
"""
from sources import *

# ==================================================
# #### Defining constants #####
# ==================================================
G = 4*np.pi**2
AU = 1.495978707E8  # (km)
M_TOT = 1
END_TIME = 100


# ====================================================
# #### Defining derivatives of the object's trajectory#####
# ====================================================
def derivatives_twobody(variables, *_, g=G, m_tot=M_TOT):
    """Calculates the derivatives with respect to time of the variables in the two body system.

    We make use of the EOM for relative motion between the two masses. Note that these EOM have no explicit time
    depencency; the second positional argument of this function is in place for compatibility with runge_kutta(),
    but we don't need any actual time input.

    Args:
        variables: NumPy array containing the relative variables x, y, vx, and vy
        g (optional): Gravitational constant G in Newton's law of gravitation.
        m_tot (optional): Total mass M1 + M2 contained in the system.

    Returns:
        A NumPy array containing the derivatives dx/dt, dy/dt, dvx/dt, and dvy/dt.
    """
    x, y, vx, vy = variables
    r = np.sqrt(x**2 + y**2)  # For use in dvx and dvy
    const = -g * m_tot / r**3   # Ditto

    # We calculate the actual derivatives
    dx = vx
    dy = vy
    dvx = const * x
    dvy = const * y

    return np.array([dx, dy, dvx, dvy])


# ==================================================
# #### Producing data file and running simulation #####
# ==================================================
def main():
    m_1 = set_value("Please enter the mass of the first primary in kg: ")
    m_2 = set_value("Please enter the mass of the second primary in kg: ")
    m1 = m_1 / (m_1 + m_2)
    m2 = m_2 / (m_1 + m_2)
    x = set_value("Please enter the initial x separation between the masses in km: ")/AU
    y = set_value("Please enter the initial y separation between the masses in km: ")/AU
    vx = set_value("Please enter the initial relative x velocity in km/s: ")*3600*24*365.2422/AU
    vy = set_value("Please enter the initial relative y velocity in km/s: ")*3600*24*365.2422/AU
    dt = set_value("Please enter the integration step dt in years (e.g. 0.01): ")
    variables = np.array([x, y, vx, vy])
    produce_data(1, variables, derivatives_twobody, END_TIME, dt, m1, m2)


main()
