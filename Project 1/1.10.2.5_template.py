#!/usr/bin/env python3
# -*- coding: utf-8 -*-


########################## PART 5 ##########################
"""
In conclusion, we can consider the case of a double pendulum, which is one of the simplest possible chaotic systems.

"""



import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
plt.close('all')

t = np.linspace(0,10,1001)

dt = t[1] - t[0] #timestep
g = 9.8 #m/s**2

#Define the parameters of the double pendulum
L1 = 1.
L2 = 1.
m1 = 1.
m2 = 1.
theta1 = np.pi/2#In radian
theta2 = np.pi#In radians


#Define the derivatives of theta and omega (the function to integrate)
def derivatives_double_pendulum(state, t):
    """
        state: numpy array giving the state of the pendulum at time t:
            theta1, omega1, theta2, omega2
        
        Returns an np.array (theta1',omega1',theta2',omega2') with the derivatives of theta and omega
    """
    #YOUR CODE GOES HERE


def runge_kutta(old_state, t, dt, derivatives):
    """
        state: numpy array giving the state of the pendulum at time t
        t: starting time
        dt: integration step
        derivatives: function that calculate the derivatives of the coordinates
        
        Function that performs an integration step using to runge-kutta algorithm
        
        Returns an np.array containing the new state (theta1, omega1, theta2, omega2)
    """
    
    #Calculate the ks
    k1 = dt * derivatives(old_state, t)
    k2 = dt * derivatives(old_state + (0.5 * k1), t + 0.5 * dt)
    k3 = dt * derivatives(old_state + (0.5 * k2), t + 0.5 * dt)
    k4 = dt * derivatives(old_state + k3, t + dt)
    
    #And consequently the new state of the system
    new_state = old_state + (k1 + 2*k2 + 2*k3 + k4) / 6.
    
    return new_state


def pendulum_location(state):
    """
    Return an array with the pendulum end locations:
    np.array([x1, y2, x2, y2])
    """
    #YOUR CODE GOES HERE


state = np.array([[theta1,0,theta2,0]]) #initial state
x_y = pendulum_location(state[0,:])

for index,time in enumerate(t[1:]): #update the pendulum angle/velocity/location using the runge kutta integration
    state_new = runge_kutta(state[index,:],time,dt,derivatives_double_pendulum)
    x_y_new = pendulum_location(state_new)
    
    state = np.vstack((state,state_new))
    x_y = np.vstack((x_y,x_y_new))
    
    
def execute_part5():
    t_window = 0.1*t[-1]
    
    fig,ax = plt.subplots(2,2)
    
    anim_theta1,      = ax[0,0].plot([],[],'b-')
    anim_theta2,      = ax[0,0].plot([],[],'r-')
    ax[0,0].set_xlim(0, t_window)
    ax[0,0].set_ylim(-1.1*np.max((np.abs(state[:,0]),np.abs(state[:,2]))),1.1*np.max((np.abs(state[:,0]),np.abs(state[:,2]))))
    ax[0,0].set_xlabel('t')
    ax[0,0].set_ylabel(r'$\theta$(t)')
    
    anim_omega1,      = ax[0,1].plot([],[],'b-')
    anim_omega2,      = ax[0,1].plot([],[],'r-')
    ax[0,1].set_xlim(0, t_window)
    ax[0,1].set_ylim(-1.1*np.max((np.abs(state[:,1]),np.abs(state[:,3]))),1.1*np.max((np.abs(state[:,1]),np.abs(state[:,3]))))
    ax[0,1].set_xlabel('t')
    ax[0,1].set_ylabel(r'$\omega$(t)')
    
    anim_pendulum,      = ax[1,0].plot([],[],'bo-')
    ax[1,0].set_xlim(-1.1*np.max(np.abs(x_y[:,2])),1.1*np.max(np.abs(x_y[:,2])))
    ax[1,0].set_ylim(-1.1*np.max(np.abs(x_y[:,3])),1.1*np.max(np.abs(x_y[:,3])))
    ax[1,0].set_xlabel('x(t)')
    ax[1,0].set_ylabel('y(t)')
    
    anim_phase1,         = ax[1,1].plot([],[],'b-')
    anim_phase2,         = ax[1,1].plot([],[],'r-')
    ax[1,1].set_xlim(-1.1*np.max((np.abs(state[:,0]),np.abs(state[:,2]))),1.1*np.max((np.abs(state[:,0]),np.abs(state[:,2]))))
    ax[1,1].set_ylim(-1.1*np.max((np.abs(state[:,1]),np.abs(state[:,3]))),1.1*np.max((np.abs(state[:,1]),np.abs(state[:,3]))))
    ax[1,1].set_xlabel(r'$\theta$(t)')
    ax[1,1].set_ylabel(r'$\omega$(t)')
    
    
    def init():
        anim_theta1.set_data([], [])
        anim_theta2.set_data([], [])
        anim_omega1.set_data([], [])
        anim_omega2.set_data([], [])
        anim_pendulum.set_data([], [])
        anim_phase1.set_data([], [])
        anim_phase2.set_data([], [])
    
    
    #The animation function which is called every frame
    def animate(i, state_,x_y_):
    
        anim_theta1.set_data(t[0:i],state_[0:i,0])
        anim_theta2.set_data(t[0:i],state_[0:i,2])
        anim_omega1.set_data(t[0:i],state_[0:i,1])
        anim_omega2.set_data(t[0:i],state_[0:i,3])
        anim_pendulum.set_data([0,x_y_[i,0],x_y_[i,2]],[0,x_y_[i,1],x_y_[i,3]])
        anim_phase1.set_data(state_[0:i,0],state_[0:i,1])
        anim_phase2.set_data(state_[0:i,2],state_[0:i,3])
        
        if t[i] > t_window:
            ax[0,0].set_xlim(t[i]-t_window, t[i])
            ax[0,1].set_xlim(t[i]-t_window, t[i])
    
    #Call the animator
    anim = animation.FuncAnimation(fig, animate, init_func=init, fargs=(state,x_y), interval=1)
    return anim    