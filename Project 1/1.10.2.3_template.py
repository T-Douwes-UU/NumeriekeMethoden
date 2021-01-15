#!/usr/bin/env python3
# -*- coding: utf-8 -*-

########################## PART 3 ##########################
"""
We can now make the pendulum a bit more realistic, by adding a viscous friction with air.
The equation of motion of the pendulum with the addition of viscous friction become
    theta''(t) = -(g/L) * sin(theta(t)) - beta * theta'(t)
with beta being a viscous friction coefficient.
The equations can be decomposed in
    omega'(t) = -(g/L) * sin(theta(t)) - beta * omega(t)
    theta'(t) = omega(t)
and can again be integrated using the runge-kutta algorithm.
The resulting trajectories are shown together with the ones of the pendulum without friction.
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
plt.close('all')

t = np.linspace(0,100,10001)

dt = t[1] - t[0] #timestep

g = 9.8 #m/s**2
L = 1.
theta_0 = np.deg2rad(60)
omega_0 = 0

#Define the viscous coefficient
beta = 0.5

def derivatives_simple_pendulum(state, t):
    """
        state: numpy array giving the state of the pendulum at time t (theta and omega)
        
        Function that returns the derivatives
            theta'(t) = omega(t)
            omega'(t) = ...
            
        Returns a np.array containing the derivatives (theta', omega')
    """

    #YOUR CODE GOES HERE

#Define the derivatives of theta and omega (the function to integrate)
def derivatives_pendulum_with_friction(state, t):
    """
        state: numpy array giving the state of the pendulum at time t
        
        Returns an np.array (theta',omega') with the derivatives of theta and omega
    """
    #YOUR CODE GOES HERE


def runge_kutta(old_state, t, dt, derivatives):
    """
        state: numpy array giving the state of the pendulum at time t
        t: starting time
        dt: integration step
        derivatives: function that calculate the derivatives of the coordinates
        
        Function that performs an integration step using to runge-kutta algorithm
        
        Returns an np.array containing the new state (theta, omega)
    """
    
    #Calculate the ks
    k1 = dt * derivatives(old_state, t)
    k2 = dt * derivatives(old_state + (0.5 * k1), t + 0.5 * dt)
    k3 = dt * derivatives(old_state + (0.5 * k2), t + 0.5 * dt)
    k4 = dt * derivatives(old_state + k3, t + dt)
    
    #And consequently the new state of the system
    new_state = old_state + (k1 + 2*k2 + 2*k3 + k4) / 6.
    
    return new_state

def pendulum_location(theta,L):
    """
    Returns the location of the pendulum end as a function of the angle and length
    """
    x = L * np.sin(theta)
    y = - L * np.cos(theta)
    return np.array([x,y])


state_fric = np.array([[theta_0,omega_0]]) #create arrays containing the initial angle/velocity/location of the pendulum
x_y_fric = pendulum_location(theta_0,L)

state_simple = np.array([[theta_0,omega_0]])
x_y_simple = pendulum_location(theta_0,L)

for index,time in enumerate(t[1:]): #update the pendulum angle/velocity/location using the runge kutta integration
    state_new_fric = runge_kutta(state_fric[index,:],time,dt,derivatives_pendulum_with_friction)
    x_y_new_fric = pendulum_location(state_new_fric[0],L)
    
    state_new_simple = runge_kutta(state_simple[index,:],time,dt,derivatives_simple_pendulum)
    x_y_new_simple = pendulum_location(state_new_simple[0],L)
    
    state_fric = np.vstack((state_fric,state_new_fric))
    x_y_fric = np.vstack((x_y_fric,x_y_new_fric))
    
    state_simple = np.vstack((state_simple,state_new_simple))
    x_y_simple = np.vstack((x_y_simple,x_y_new_simple))


def execute_part3():
    t_window = 0.1*t[-1]
    
    fig,ax = plt.subplots(2,2)
    
    anim_theta_fric,         = ax[0,0].plot([],[],'b-')
    anim_theta_simple,      = ax[0,0].plot([],[],'b--')
    ax[0,0].set_xlim(0, t_window)
    ax[0,0].set_ylim(-1.1*np.max(np.abs(state_simple[:,0])),1.1*np.max(np.abs(state_simple[:,0])))
    ax[0,0].set_xlabel('t')
    ax[0,0].set_ylabel(r'$\theta$(t)')
    
    anim_omega_fric,         = ax[0,1].plot([],[],'r-')
    anim_omega_simple,      = ax[0,1].plot([],[],'r--')
    ax[0,1].set_xlim(0, t_window)
    ax[0,1].set_ylim(-1.1*np.max(np.abs(state_simple[:,1])),1.1*np.max(np.abs(state_simple[:,1])))
    ax[0,1].set_xlabel('t')
    ax[0,1].set_ylabel(r'$\omega$(t)')
    
    anim_pendulum,      = ax[1,0].plot([],[],'bo-')
    anim_trajectory,    = ax[1,0].plot([],[],'r-')
    ax[1,0].set_xlim(-1.1*np.max(np.abs(x_y_simple[:,0])),1.1*np.max(np.abs(x_y_simple[:,0])))
    ax[1,0].set_ylim(-1.1*np.max(np.abs(x_y_simple[:,1])),1.1*np.max(np.abs(x_y_simple[:,1])))
    ax[1,0].set_xlabel('x(t)')
    ax[1,0].set_ylabel('y(t)')
    
    anim_phase,         = ax[1,1].plot([],[],'b-')
    ax[1,1].set_xlim(-1.1*np.max(np.abs(state_simple[:,0])),1.1*np.max(np.abs(state_simple[:,0])))
    ax[1,1].set_ylim(-1.1*np.max(np.abs(state_simple[:,1])),1.1*np.max(np.abs(state_simple[:,1])))
    ax[1,1].set_xlabel(r'$\theta$(t)')
    ax[1,1].set_ylabel(r'$\omega$(t)')
    
    
    
    def init():
        anim_theta_fric.set_data([], [])
        anim_theta_simple.set_data([], [])
        anim_omega_fric.set_data([], [])
        anim_omega_simple.set_data([], [])
        anim_pendulum.set_data([], [])
        anim_trajectory.set_data([], [])
        anim_phase.set_data([], [])
        
    
    #The animation function which is called every frame
    def animate(i, state,state_sa,x_y,x_y_sa):
    
        anim_theta_fric.set_data(t[0:i],state_fric[0:i,0])
        anim_theta_simple.set_data(t[0:i],state_simple[0:i,0])
        anim_omega_fric.set_data(t[0:i],state_fric[0:i,1])
        anim_omega_simple.set_data(t[0:i],state_simple[0:i,1])
        anim_pendulum.set_data([0,x_y_fric[i,0]],[0,x_y_fric[i,1]])
        anim_trajectory.set_data(x_y_fric[0:i,0],x_y_fric[0:i,1])
        anim_phase.set_data(state_fric[0:i,0],state_fric[0:i,1])
        
        if t[i] > t_window:
            ax[0,0].set_xlim(t[i]-t_window, t[i])
            ax[0,1].set_xlim(t[i]-t_window, t[i])
    
    #Call the animator
    anim = animation.FuncAnimation(fig, animate, init_func=init, fargs=(state_fric,state_simple,x_y_fric,x_y_simple), interval=10) #change interval for the plotting speed
    return anim    
    

