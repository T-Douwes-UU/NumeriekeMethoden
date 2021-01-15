"""
This file produces the data for Task 3. Running it prompts the user to input initial values for the system,
after which the necessary data is produced and output to a file in the parent directory.

This file must be in the same folder as sources.py to be run successfully.
"""
from sources import *

# =============================
# #### Defining constants #####
# =============================
MU_E = 2  # Mean number of nucleons per electron
M_E = 9.10938370e-31  # Electron mass in kg
M_H = 1.6735575e-27  # Hydrogen mass in kg
C = 299792458.  # Speed of light in m/s
H = 6.62607004e-34  # Planck's constant in m^2*kg/s
G = 6.67430e-11  # Gravitational constant in m^3/kg*s^2

A = 8/3*np.pi * (M_E*C/H)**3 * M_H
P_0 = 8/3*np.pi * (M_E*C/H)**3 * M_E*C**2
ALPHA_R = np.sqrt(P_0 / (4*np.pi*G)) / (MU_E*A)
ALPHA_M = np.sqrt((P_0 / G)**3 / (4*np.pi)) / (MU_E*A)**2


# ======================================================
# #### Calculating starting values for integration #####
# ======================================================
def starting_values(rho_c, dxi):
    """Produces the starting values x(dxi) and mu(dxi) for numerical integration of the structure of a white dwarf.

    Args:
        rho_c: The core density of the white dwarf in kg/m^3.
        dxi: Integration step for xi (dimensionless).

    Returns:
        A NumPy array containing x(dxi) and mu(dxi).
    """
    x_c = (rho_c / (A*MU_E))**(1./3)

    x_dxi = x_c - (x_c**2)/6 * np.sqrt(1 + x_c**2) * dxi**2
    mu_dxi = (x_c**3)/3 * dxi**3
    return np.array([x_dxi, mu_dxi])


# ===============================
# #### Defining derivatives #####
# ===============================
def derivatives_whitedwarf(variables, xi):
    """Calculates the derivatives with respect to xi of the variables in the structure equations of a white dwarf star.

    We make use of the equations using dimensionless variables x and mu. We integrate over xi, not time, but we can
    use our Runge Kutta implementation just fine regardless.

    Args:
        variables: NumPy array containing the variables x and mu
        xi: Dimensionless variable representing the distance from the center of the white dwarf.

    Returns:
        A NumPy array containing the derivatives dx/dxi and dmu/dxi.
    """
    x, mu = variables

    dx = -(mu / xi**2) * (np.sqrt(1 + x**2) / x)
    dmu = x**3 * xi**2

    return np.array([dx, dmu])


# ======================================================
# #### Producing data file and running integration #####
# ======================================================
def produce_data_task3(rho_c, dxi, end_condition=1e-4, maxlines=10_000):
    """
    Produces the file containing the simulated data for this task. Each line contains a value of xi, x and mu.

    Args:
        rho_c: The core density of the white dwarf in kg/m^3.
        dxi: Integration step for xi (dimensionless).
        end_condition: Integration ends once x drops below this value.
        maxlines: Maximum amount of lines to write out. This is more of a measure to prevent the function from running
            forever in case of a mistake or improper parameters.

    Returns:
        A text file named 'data_task3_<date>_<time>.txt'.
    """
    print("Working...")

    xi = dxi
    variables = starting_values(rho_c, dxi)

    now = datetime.now()  # Get the current date and time
    filename = now.strftime(f"data_task3_%Y-%m-%d_%H.%M.%S.txt")
    with open(filename, 'w') as f:
        f.write(f"Task 3 data, written out at {now.strftime('%H:%M:%S on %d-%m-%Y')}.\n"
                f"Core density: {rho_c} kg/m^3. Integration step dxi: {dxi}\n\n"
                f"Data starts on line 6 and is formatted as <xi (tab) x (tab) mu>.\n\n")

        lines = 5
        halves = 0

        while variables[0] > end_condition and lines < maxlines:
            [x, mu] = variables
            f.write(f"{xi}\t{x}\t{mu}\n")
            new_vars = runge_kutta(variables, derivatives_whitedwarf, xi, dxi)

            # Handling two unwanted scenarios by increasing precision:
            # - The Runge-Kutta method calculates ahead beyond x<0, causing dx>0
            # - The integration step takes us beyond x<0 in one step (we want 0<x<end_condition)
            while new_vars[0] > x or new_vars[0] < 0:
                dxi /= 2  # Integration step is halved
                halves += 1  # Bookkeeping
                new_vars = runge_kutta(variables, derivatives_whitedwarf, xi, dxi)

            variables = new_vars
            xi += dxi
            lines += 1

        [x, mu] = variables
        f.write(f"{xi}\t{x}\t{mu}\n")  # To include the last data point
    if lines < maxlines:
        print(f"Finished writing data to {filename}. Halved integration step {halves} times.")
    else:
        print(f"Maximum line number reached, stopped writing to {filename}.")


def main():
    rho_c = 10**set_value("Please enter the order of the core pressure (base 10 log of rho_c in kg/m^3): ")
    dxi = set_value("Please enter the integration step of r in km: ") * 1000 / ALPHA_R
    produce_data_task3(rho_c, dxi, end_condition=1e-4, maxlines=10_000)


main()
