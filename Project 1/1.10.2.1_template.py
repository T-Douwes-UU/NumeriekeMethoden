#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
In these exercises, we want to apply the Runge-Kutta numerical integration method to simple but interesting problems,
of increasing complexity.
We will study the behaviour of pendulums, starting from the easiest possible case: a simple pendulum, with no friction,
in the small angles approximation. For this case, an analytical solution exists which we can plot and later use as a test.
We will then explore the condition of bigger angles, for which an analytical solution does not exist, numerically
integrating the equations of motion.
We will then introduce the viscous friction of the pendulum with air, dumping the equations of motion, and a forcing
pulsation.
In the end, we will consider the case of coupled pendulums, which is an extremely simple but appealing chaotic system.
"""


########################## PART 1 ##########################
"""
We can start by considering the simple pendulum of length L.
The status of a pendulum can be entirely described by its angle theta with the vertical direction and its angular 
velocity omega.
The equation of motion is the following
    theta''(t) = -(g/L) * sin(theta(t))
    
In the approximation of small angles, such that sin(theta) can be replaced by theta, an analytical solution exists.
Derive this solution, and plug it into the small_angle_approximation function
    
We can then plot the evolution of theta (in blue) and omega (in red) with time for a certain starting condition
theta0.
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
plt.close('all')


t = np.linspace(0,10,1001)

g = 9.8 #m/s**2
L = 1.
theta_0 = np.deg2rad(60)
omega_0 = 0

#Function that returns the theta and the omega of the simple pendulum at time t in the small angle approximation,
#given the starting angle and the time
def small_angle_approximation(t, theta0):
    """
        t: time
        theta0: starting angle
        
    Function that returns the angle and the velocity of a simple pendulum in the small angle approximation.
    
    Returns the tuple (theta(t), omega(t)) containing the angle and the angular velocity of the pendulum
    """
    f = np.sqrt(g/L)
    thetat = theta0 * np.cos(f*t)
    omegat = -theta0 * f * np.sin(f*t)
    return thetat, omegat


def pendulum_location(theta,L):
    """
    Returns the location of the pendulum end as a function of the angle and length
    """
    x = L * np.sin(theta)
    y = - L * np.cos(theta)
    return np.array([x,y])


theta = np.array([theta_0]) #create arrays containing the initial angle/velocity/location of the pendulum
omega = np.array([omega_0])
x_y = pendulum_location(theta_0,L)


for time in t[1:]: #update the pendulum angle/velocity/location in time using the small angle approximation
    # question: why do we use t[1:] here?
    
    theta_new,omega_new = small_angle_approximation(time,theta_0)
    x_y_new = pendulum_location(theta_new,L)
    
    theta = np.append(theta,theta_new) #append the new values to the arrays
    omega = np.append(omega,omega_new)
    x_y = np.vstack((x_y,x_y_new)) #try typing 'help(np.vstack)' in your console if you wonder what this line here is doing

    
def execute_part1():
    """
    Animate your results (type 'execute_part1()' in your iPython console)
    
    Try to get a rough idea what is going on here if you have some time left! Understanding this 
    will come in very handy in the second part of numerical methods in 4 weeks, where we don't provide
    you anymore with these nice templates!
    """
    t_window = 0.1*t[-1] #you can play around with this variable if you want to plot e.g. a longer time window
    
    fig,ax = plt.subplots(2,2)
    
    anim_theta,         = ax[0,0].plot([],[],'b-')
    ax[0,0].set_xlim(0, t_window)
    ax[0,0].set_ylim(-1.1*np.max(np.abs(theta)),1.1*np.max(np.abs(theta)))
    ax[0,0].set_xlabel('t')
    ax[0,0].set_ylabel(r'$\theta$(t)')
    
    anim_omega,         = ax[0,1].plot([],[],'r-')
    ax[0,1].set_xlim(0, t_window)
    ax[0,1].set_ylim(-1.1*np.max(np.abs(omega)),1.1*np.max(np.abs(omega)))
    ax[0,1].set_xlabel('t')
    ax[0,1].set_ylabel(r'$\omega$(t)')
    
    anim_pendulum,      = ax[1,0].plot([],[],'bo-')
    anim_trajectory,    = ax[1,0].plot([],[],'r-')
    ax[1,0].set_xlim(-1.1*np.max(np.abs(x_y[:,0])),1.1*np.max(np.abs(x_y[:,0])))
    ax[1,0].set_ylim(-1.1*np.max(np.abs(x_y[:,1])),1.1*np.max(np.abs(x_y[:,1])))
    ax[1,0].set_xlabel('x(t)')
    ax[1,0].set_ylabel('y(t)')
    
    anim_phase,         = ax[1,1].plot([],[],'b-')
    ax[1,1].set_xlim(-1.1*np.max(np.abs(theta)),1.1*np.max(np.abs(theta)))
    ax[1,1].set_ylim(-1.1*np.max(np.abs(omega)),1.1*np.max(np.abs(omega)))
    ax[1,1].set_xlabel(r'$\theta$(t)')
    ax[1,1].set_ylabel(r'$\omega$(t)')
    
    
    
    def init():
        anim_theta.set_data([], [])
        anim_omega.set_data([], [])
        anim_pendulum.set_data([], [])
        anim_trajectory.set_data([], [])
        anim_phase.set_data([], [])
        
    
    #The animation function which is called every frame
    def animate(i, theta,omega,x_y):
    
        anim_theta.set_data(t[0:i],theta[0:i])
        anim_omega.set_data(t[0:i],omega[0:i])
        anim_pendulum.set_data([0,x_y[i,0]],[0,x_y[i,1]])
        anim_trajectory.set_data(x_y[0:i,0],x_y[0:i,1])
        anim_phase.set_data(theta[0:i],omega[0:i])
        
        if t[i] > t_window:
            ax[0,0].set_xlim(t[i]-t_window, t[i])
            ax[0,1].set_xlim(t[i]-t_window, t[i])
    
    #Call the animator
    anim = animation.FuncAnimation(fig, animate, init_func=init, fargs=(theta,omega,x_y), interval=2) #change interval for the plotting speed
    return anim