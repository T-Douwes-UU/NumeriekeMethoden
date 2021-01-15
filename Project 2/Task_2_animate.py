"""
This file loads the data from the Task_2_data.py script for Task 2 and produces
an animation.

This file must be in the same folder as sources.py to be run successfully.
"""
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from sources import *

# ==================================================
# #### Defining constants #####
# ==================================================
G = 1
omega = 1


# ==================================================
# #### Animating data #####
# ==================================================
def animate_data_three_body():
    """
    Takes the data from the data file and animates it to give a visual representation
    of the simulation. Also outputs plots of the total energy of the asteroid.
    """
    # Obtaining data and calculating the energy
    t, x, y, m1, m2, vx, vy = load_data()
    x1, x2, y1, y2 = positions(1, 0, m1, m2)
    R1 = np.sqrt((x - x1)**2 + (y - y1)**2)
    R2 = np.sqrt((x - x2)**2 + (y - y2)**2)
    v = np.sqrt(vx**2 + vy**2)
    E1 = 0.5 * v**2 - G * (m1 / R1 + m2 / R2)  # Kinetic and Potential energy / mass asteroid
    E_rot = - 0.5 * omega**2 * (x**2 + y**2)  # Rotational energy / mass asteroid
    E2 = E1 + E_rot  # Energy with rotational energy included
    print(f"Found average energy is: E = {np.average(E2)}")

    # Plotting found energy values without rotational energy
    fig1, ax_E1 = plt.subplots(1)
    ax_E1.set_xlim(0, max(t))
    ax_E1.set_xlabel(r'$t$ (yr)')
    ax_E1.set_ylabel(r'$E/m$')
    ax_E1.set_title('Energy per mass of Trojan asteroid in horseshoe orbit')
    ax_E1.plot(t, E1, color='olivedrab', label='Trojan asteroid')
    fig1.savefig("Three-body Energy.png", dpi=300, bbox_inches='tight')

    # Plotting found energy values with the rotational energy included
    fig2, ax_E2 = plt.subplots(1)
    ax_E2.set_xlim(0, max(t))
    ax_E2.set_ylim(-1.51, -1.49)
    ax_E2.set_xlabel(r'$t$ (yr)')
    ax_E2.set_ylabel(r'$E/m$')
    ax_E2.set_title('Energy per mass of Trojan asteroid in horseshoe orbit (incl. $E_{rot}$)')
    ax_E2.plot(t, E2, color='olivedrab', label='Trojan asteroid')
    fig2.savefig("Three-body Energy (incl Erot).png", dpi=300, bbox_inches='tight')

    # Setting up the figure for the animation and plotting the stationary parts
    fig3, ax = plt.subplots(1)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Three body System')
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    # plotting the trajectory of the small object
    ax.plot(x, y, color='cornflowerblue', lw=0.5, label='Trajectory Trojan asteroid')
    ax.plot(x1, y1, color='goldenrod', linestyle='None', marker='o',
            markersize=10, label='First primary body (Sun)')
    ax.plot(x2, y2, color='seagreen', linestyle='None', marker='o',
            markersize=6, label='Second primary body (Jupiter)')
    ax.plot(-(m1-m2)/2, np.sqrt(3)/2, color='red', linestyle='None',
            marker='o', markersize=2, label='L4 Lagrange point')
    ax.plot(-(m1-m2)/2, -np.sqrt(3)/2, color='chocolate', linestyle='None',
            marker='o', markersize=2, label='L5 Lagrange point')
    ax.legend(fontsize='xx-small')
    fig3.savefig("Three-body problem.png", dpi=300, bbox_inches='tight')
    lasteroid, = ax.plot([], [], 'o', color='magenta', marker='o',
                         markersize=4, label='Trojan asteroid')
    plt.legend()

    # Setting up the data for both objects in the animation, fixing the big object
    def animate(i):
        lasteroid.set_data(x[i], y[i])
        return lasteroid,

    return animation.FuncAnimation(fig3, animate, frames=len(x), interval=40, blit=True)


anim = animate_data_three_body()
