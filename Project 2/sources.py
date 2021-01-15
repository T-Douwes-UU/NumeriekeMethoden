"""
This file stores functions that are used across multiple tasks, in order of usage.
"""
from datetime import datetime
import numpy as np
import re
from tkinter import Tk
from tkinter.filedialog import askopenfilename

Tk().withdraw()


def set_value(prompt):
    """Simple input function with error handling for non-number entries.

    Args:
        prompt: A string containing an input prompt.

    Returns:
        The number that is put in, as a float.
    """
    while True:
        try:
            value = float(input(prompt))
            return value
        except ValueError:
            print("Not a number, please try again.")


def runge_kutta(old_state, derivatives, t, dt):
    """Performs an integration step using the Runge-Kutta algorithm.

    Note that this method is defined in terms of a time variable t, but it works just as well for other
    types of variables.

    Args:
        old_state: NumPy array giving the state of the system variables at time t
        derivatives: function that calculates the derivatives of the coordinates
        t: starting time
        dt: integration step

    Returns:
        A NumPy array containing the new state of the system variables at time t+dt
    """
    # We calculate the ks
    k1 = dt * derivatives(old_state, t)
    k2 = dt * derivatives(old_state + (0.5 * k1), t + 0.5 * dt)
    k3 = dt * derivatives(old_state + (0.5 * k2), t + 0.5 * dt)
    k4 = dt * derivatives(old_state + k3, t + dt)

    # And consequently the new state of the system
    new_state = old_state + (k1 + 2.*k2 + 2.*k3 + k4) / 6.

    return new_state


def produce_data(tasknumber, variables, derivatives, end_time, dt, m1, m2):
    """
    Produces the file containing the simulated data for tasks one and two. Each line contains the time,
    (relative) positions and (relative) velocities.

    Args:
        tasknumber: Number of the task to produce data for.
        variables: NumPy array containing the initial variables x, y, vx and vy.
        derivatives: Function that calculates the derivatives for the specific task.
        end_time: The time at which the simulation ends.
        dt: Time step used for numerical integration.
        m1: The large mass.
        m2: The small mass.

    Returns:
        A text file named 'data_task<tasknumber>_<date>_<time>.txt'.
    """
    time = np.arange(0, end_time + dt, dt)
    print("Working...")

    now = datetime.now()  # Get the current date and time
    filename = now.strftime(f"data_task{tasknumber}_%Y-%m-%d_%H.%M.%S.txt")
    with open(filename, 'w') as f:
        f.write(f"Task {tasknumber} data, written out at {now.strftime('%H:%M:%S on %d-%m-%Y')}.\n"
                f"Masses: {m1}, {m2}\n"
                f"Time runs from t=0 to t={end_time} in {len(time)} steps of dt={dt}\n"
                f"Data starts on line 6 and is formatted as <t (tab) x (tab) y (tab) vx (tab) vy>.\n\n")

        for t in time:
            [x, y, vx, vy] = variables
            f.write(f"{t}\t{x}\t{y}\t{vx}\t{vy}\n")
            variables = runge_kutta(variables, derivatives, t, dt)
    print(f"Finished writing data to {filename}.")


def positions(x, y, m1, m2):
    """
    Takes the seperation between the two masses and calculates the positions of both
    masses.
    
    Args:
        x: Seperation between the masses in the x direction.
        y: Seperation between the masses in the y direction.
        m1: The fraction of the total mass made up by m_1. (m_1 / m_tot)
        m2: The fraction of the total mass made up by m_2. (m_2 / m_tot)
    
    Returns:
        x1, x2, y1, y2
    """
    x1 = x * m2
    x2 = -x * m1
    y1 = -y * m2
    y2 = y * m1
    return x1, x2, y1, y2


def load_data():
    """
    Loads in the file containing the simulated data. unpacking it into the time, relative position and
    relative velocities.

    Returns:
        A tuple containing: arrays for the time and positions; the two masses as floats; arrays for the velocities.
    """
    print("Please choose a data file to animate.")
    filename = askopenfilename(title="Choose a file", filetypes=[("Text files", ".txt")])
    with open(filename) as f:
        m1, m2 = re.split(': |, ', f.readlines()[1])[-2:]
    data = np.genfromtxt(filename, skip_header=4, delimiter='\t')
    t = data[:, 0]
    x = data[:, 1]
    y = data[:, 2]
    vx = data[:, 3]
    vy = data[:, 4]
    return t, x, y, float(m1), float(m2), vx, vy
