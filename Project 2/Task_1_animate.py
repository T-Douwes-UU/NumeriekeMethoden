"""
This file loads the data from the Task_1_data.py script for Task 1 and produces
an animation.

This file must be in the same folder as sources.py to be run successfully.
"""
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from sources import *

# ==================================================
# #### Defining constants #####
# ==================================================
G = 4*np.pi**2


# ==================================================
# #### Animating data #####
# ==================================================
def animate_data_two_body():
    """
    Takes the data from the data file and animates it to give a visual representation
    of the simulation. In addition, it also calculates and prints out the period of the orbit.
    """
    # obtaining data
    t, x, y, m1, m2 = load_data()[0:5]
    x1, x2, y1, y2 = positions(x, y, m1, m2)

    # calculating period
    y_old = 0
    period_pos = np.array([])
    for i in range(len(y)):
        if y_old < 0 <= y[i]:
            period_pos = np.append(period_pos, i)
        y_old = y[i]
    period = t[int(period_pos[-1])] / len(period_pos)
    period_kepler = np.sqrt(4*np.pi**2*abs(x[0])**3/(G*(m1+m2)))
    print(f"The period of this orbit calculated out of data is P = {period} years")
    print(f"The period according to Keplers law is P = {period_kepler} years")

    # setting up the figure
    fig, ax = plt.subplots(1)
    ax.set_xlim(-1.5*abs(x[0]), 1.5*abs(x[0]))
    ax.set_ylim(-1.5*abs(x[0]), 1.5*abs(x[0]))

    # plotting the trajectory of the small object
    plt.plot(x2, y2, color='lightskyblue', lw=0.5, label='trajectory secondary body')
    # setting up the plots of the bodies
    lprimary, = ax.plot([], [], 'o', color='red', marker='o',
                        markersize=14, label='primary body')
    lsecondary, = ax.plot([], [], 'o', color='blue', marker='o',
                          markersize=10, label='secondary body')
    # Report colour guide:
    # - Sun: color='goldenrod'
    # - Earth: color='royalblue'; trajectory: color='lightskyblue'
    # - Mercury: color='gray'; trajectory: color='gainsboro'
    # - Jupiter: color='chocolate'; trajectory: color='lightsalmon'
    plt.xlabel('x (AU)')
    plt.ylabel('y (AU)')
    plt.title('Two-body system')
    plt.legend(fontsize='x-small')

    # initiating the datasets for both objects in the animation
    def init():
        lprimary.set_data([], [])
        lsecondary.set_data([], [])
        return lprimary, lsecondary
    
    # setting up the data for both objects in the animation, fixing the big object
    def animate(frame):
        lprimary.set_data(x1[frame], y1[frame])
        lsecondary.set_data(x2[frame], y2[frame])
        return lprimary, lsecondary

    # Uncomment these lines to export a static .png as in the report
    # plt.plot(x1[0], y1[0], 'o', color='red', marker='o', markersize=14, label='primary body (the Sun)')
    # plt.plot(x2[0], y2[0], 'o', color='blue', marker='o', markersize=10, label='secondary body')
    # plt.savefig('Two-body system.png', dpi=300, bbox_inches='tight')
    
    return animation.FuncAnimation(fig, animate, frames=len(x), init_func=init, interval=2)


anim = animate_data_two_body()
