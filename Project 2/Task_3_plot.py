"""
This file loads the data from the Task_3_data.py script for Task 3 and produces
plots of the data.

This file must be in the same folder as sources.py to be run successfully.
"""
import matplotlib.pyplot as plt
from sources import *

Tk().withdraw()

# ==================================================
# #### Defining constants #####
# ==================================================
MU_E = 2  # Mean number of nucleons per electron
M_SUN = 1.9885e30  # Solar mass in kg
M_E = 9.10938370e-31  # Electron mass in kg
M_H = 1.6735575e-27  # Hydrogen mass in kg
C = 299792458.  # Speed of light in m/s
H = 6.62607004e-34  # Planck's constant in m^2*kg/s
G = 6.67430e-11  # Gravitational constant in m^3/kg*s^2

A = 8/3*np.pi * (M_E*C/H)**3 * M_H
P_0 = 8/3*np.pi * (M_E*C/H)**3 * M_E*C**2
ALPHA_R = np.sqrt(P_0 / (4*np.pi*G)) / (MU_E*A)
ALPHA_M = np.sqrt((P_0 / G)**3 / (4*np.pi)) / (MU_E*A)**2


# ==================================================
# #### Plotting data #####
# ==================================================
def plot_data():
    """
    plots the density against the radius of the white dwarf and the maximum radius
    against the corresponding mass of the white dwarf using the data from the selected
    datafile.
    
    Returns: 
        a tuple containing the maximal r, rho and M value.
    """
    print("Please choose all data files you want to use.")
    filenames = askopenfilename(title="Choose all files you wish to plot", 
                                filetypes=[("Text files", ".txt")], multiple=True)
    print("Working...")
    
    # Setting up the figure for the density against radius plot:
    fig1, ax1 = plt.subplots(1)
    max_x1 = 0
    max_y1 = 0
    min_y1 = 1
    
    # Setting up the figure for max radius against Mass plot:
    fig2, ax2 = plt.subplots(1)
    max_x2 = 0
    max_y2 = 0
    
    # Evaluating every data file:
    for filename in filenames:
        # Retrieving data and values for file i
        data = np.genfromtxt(filename, skip_header=4, delimiter='\t')
        with open(filename) as file:
            rho_c, dxi = [float(s) for s in re.findall(r'\d*\.\d+', file.readlines()[1])]
        xi = data[:, 0]
        x = data[:, 1]
        mu = data[:, 2]
        
        # Calculating appropriate values:
        rho = A * MU_E * x**3
        r = ALPHA_R * xi / 1000
        M = ALPHA_M * mu / M_SUN  # mass of white dwarf in Solar mass
        
        # Plotting the density plot of file i
        ax1.plot(r, rho, label=rf"$\rho_c$={rho_c:.0E}, $dr$={dxi * ALPHA_R / 1000:.0f}")
        
        if max(r) > max_x1:
            max_x1 = max(r)
        if max(rho) > max_y1:
            max_y1 = max(rho)
        if min(rho) < min_y1:
            min_y1 = min(rho)
        
        # Plotting the radius plot of file i
        ax2.plot(M[len(M)-1], r[len(r)-1], marker='o', linestyle='none',
                 label=rf"$\rho_c$={rho_c:.0E}, $dr$={dxi * ALPHA_R / 1000:.0f}")
        
        if max(M) > max_x2:
            max_x2 = max(M)
        if max(r) > max_y2:
            max_y2 = max(r)
        
    # formatting the figure for the density against radius plot:
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(right=1.2*max_x1)
    ax1.set_ylim(0.5*min_y1, 3*max_y1)
    ax1.set_xlabel(r'r ($km$)', fontsize=10)
    ax1.set_ylabel(r'$\rho$ ($kg/m^3$)', fontsize=10)
    ax1.set_title('White dwarf density plotted against radius', fontsize=12)
    ax1.legend(fontsize='xx-small')
    fig1.savefig("Task 3 density plot.png", dpi=300, bbox_inches='tight')
    
    # formatting the figure for the maximal r against corresponding M plot:
    ax2.set_xlim(0, 1.01*max_x2)
    ax2.set_ylim(0, 1.02*max_y2)
    ax2.set_xlabel(r"m ($M_\odot$)", fontsize=10)
    ax2.set_ylabel(r'r ($km$)', fontsize=10)
    ax2.set_title('White dwarf radius plotted against mass', fontsize=12)
    ax2.legend(fontsize='xx-small')
    fig2.savefig("Task 3 radius plot.png", dpi=300, bbox_inches='tight')
    
    print("Finished plotting all data.")


plot_data()
