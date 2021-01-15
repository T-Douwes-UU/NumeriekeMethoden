"""
We intend to apply Monte Carlo methods to the evaluation of integrals.

Monte Carlo methods can be directly employed to the evaluation of integrals of known functions, and therefore
to the calculation of averages with respect to fully known distributions.


We start with the evaluation of the integral
    I = int_0^1 dx exp(-x)
We will calculate this integral with deterministic numerical methods, as the methods of rectangles, trapezoids, and 
the Simpson method. Then, we will calculate it using a Monte Carlo method, both the simple one and applying importance
sampling.
"""


"""
The integral
    I = int_0^1 dx exp(-x)
can be numerically calculated employing the simple rectangles method.
In the simple rectangles method, to calculate an integral
    I = int_a^b dx f(x)
a mesh of the interval [a,b] of N+1 points xi (with x(0) = a and x(N) = b) is considered, which slices the interval
in N different sectors. Then, the area below the function f(x) in each of these slices [x(i),x(i+1)) is evaluated
by approximated the function f(x) as a constant within the slice: f(x) = f(x(i)), x in [x(i),x(i+1)).
Namely, the integral is approximated as
    I = sum_i=0^N ((b-a)/N) * f(x(i))
"""

import numpy as np
import matplotlib.pyplot as plt




#Method that performs the standard rectangle integration
def rectangles_method_integration(integrand, interval, N):
    """
        integrand: integrand function
        interval: tuple containing the extremes (a,b) of the integration interval
        N: number of slices for the integration
        
    Routine that performs the numerical integration of an integrand function in a certain interval (a,b) with the
    simple method of the rectangle.
    """
    
    #YOUR CODE GOES HERE


def execute_part1():
    #Define the parameters of the integration mesh
    N = #YOUR CODE GOES HERE
    a = #YOUR CODE GOES HERE
    b = #YOUR CODE GOES HERE
    
    #Define and the integrand function
    def integrand(x):
        #YOUR CODE GOES HERE
    
    #Perform the integration
    I = rectangles_method_integration(integrand, (a,b), N)
    trueI = 0.63212
    print('The rectangle method returns I = {0:1.6f}.\nThe calculation yields a relative error of {1:1.6f}%.\n'.format(I, (100. * ((I - trueI) / trueI))))


"""
The same integral I can be numerically calculated employing the trapezoids method.
In the trapezoid method, to calculate an integral
    I = int_a^b dx f(x)
a mesh of the interval [a,b] of N+1 points xi (with x(0) = a and x(N) = b) is considered, which slices the interval
in N different sectors. Then, the area below the function f(x) in each of these slices [x(i),x(i+1)) is evaluated
by approximated the function f(x) as a line between the points (x(i),f(x(i)) and (x(i+1),f(x(i+1))
Namely, the integral is approximated as
    I = sum_i=0^N ((b-a)/N) * ((f(x(i)) + f(x(i+1))) / 2) 
"""

#Method that performs the standard rectangle integration
def trapezoids_method_integration(integrand, interval, N):
    """
        integrand: integrand function
        interval: tuple containing the extremes (a,b) of the integration interval
        N: number of slices for the integration
        
    Routine that performs the numerical integration of an integrand function in a certain interval (a,b) with the
    trapezoids method.
    """
    
    #YOUR CODE GOES HERE


def execute_part2():
    #Define the parameters of the integration mesh
    N = #YOUR CODE GOES HERE
    a = #YOUR CODE GOES HERE
    b = #YOUR CODE GOES HERE
    
    #Define and the integrand function
    def integrand(x):
        #YOUR CODE GOES HERE
    
    #Perform the integration
    I = trapezoids_method_integration(integrand, (a,b), N)
    trueI = 0.63212
    print('The trapezoids method returns I = {0:1.6f}.\nThe calculation yields a relative error of {1:1.6f}%.\n'.format(I, (100. * ((I - trueI) / trueI))))


"""
Also, the Simpson method can be applied to the numerical evaluation of the integral.
In the Simpson method, to calculate an integral
    I = int_a^b dx f(x)
a mesh of the interval [a,b] of N+1 points xi (with x(0) = a and x(N) = b) is considered. 
Then, for each set of three points xi-1, xi, xi+1 with i odd the function is approximated as a parabola passing
through the points, and the area below it consequently evaluated. The integral is therefore approximated as
    I = sum_{i=0, i odd}^N-1 ((f(x(i-1)) + 4 * f(x(i)) f(x(i+1))) * ((b-a) / (3*N)))
"""




#Method that performs the standard rectangle integration
def simpson_method_integration(integrand, interval, N):
    """
        integrand: integrand function
        interval: tuple containing the extremes (a,b) of the integration interval
        N: number of slices for the integration
        
    Routine that performs the numerical integration of an integrand function in a certain interval (a,b) with the
    Simpson method.
    """
    
    #YOUR CODE GOES HERE
    

def execute_part3():
    #Define the parameters of the integration mesh
    N = #YOUR CODE GOES HERE
    a = #YOUR CODE GOES HERE
    b = #YOUR CODE GOES HERE
    
    assert (N%2 == 0)
    
    #Define the integrand function
    def integrand(x):
        #YOUR CODE GOES HERE
    
    I = simpson_method_integration(integrand, (a,b), N)
    trueI = 0.63212
    print('The Simpson method returns I = {0:1.6f}.\nThe calculation yields a relative error of {1:1.6f}%.\n'.format(I, (100. * ((I - trueI) / trueI))))


"""
In the end, we can employ a Monte Carlo method to numerically estimate the integral.
"""

import random
from matplotlib import animation
import collections




#Method that performs the standard rectangle integration
def MC_method_integration(integrand, interval, N, samples, squaredsamples):
    """
        integrand: integrand function
        interval: tuple containing the extremes (a,b) of the integration interval
        N: number of random numbers generated
        samples, squaredsamples: current available samples and squares of the samples
        
    Routine that performs the numerical integration of an integrand function in a certain interval (a,b) with the
    Monte Carlo method.
    Return the tuple (I, sigma, samplesn, squaredsamplesn), where:
        - I is the current value of the integral
        - sigma is the current integration error
        - samplesn and squaredsamplesn are the updated list of samples and their squares
    """
    
    #YOUR CODE GOES HERE



def execute_part4():
    
    
    Npercycle = #YOUR CODE GOES HERE #amount of MC samples per cycle
    Ncycles = #YOUR CODE GOES HERE #total amount of cycles
    a = #YOUR CODE GOES HERE
    b = #YOUR CODE GOES HERE
    
    #Define the integrand function
    def integrand(x):
        #YOUR CODE GOES HERE
    
    
    # initialize the arrays which will store the results:
    n = np.array([]) #iterations
    I = np.array([]) #integrand
    sigma = np.array([]) #integration error
    
    samples = np.array([]) #array storing the samples
    squaredsamples = np.array([])
    
    for i in range(Ncycles):
        
        (In, sigman, samplesn, squaredsamplesn) = MC_method_integration(integrand, (a,b), Npercycle, samples, squaredsamples)
    
        n = np.append(n,i)
        I = np.append(I,In)
        sigma = np.append(sigma,sigman)
    
        samples = samplesn
        squaredsamples = squaredsamplesn
 
    
    #First set up the figure, the axis, and the plot element we want to animate
    stdsize = plt.rcParams.get('figure.figsize')
    fig = plt.figure(figsize=(2*stdsize[0],stdsize[1]))
    axI = fig.add_subplot(121)
    axI.set_xlim(0, Ncycles)
    axI.set_ylim(0.60,0.65)
    lI, = axI.plot([],[], 'b-')
    axsigma = fig.add_subplot(122)
    axsigma.set_xlim(0, Ncycles)
    axsigma.set_ylim(0.,1.)
    lsigma, = axsigma.plot([],[],'r-')
    # YOU CAN ADJUST THE LABELS / LIMITS OF THE FIGURE YOURSELF
    
    
    def init():
        lI.set_data([], [])
        lsigma.set_data([],[])
    
    #And the animation function which is called every frame
    def animate(i, I, sigma):

        lI.set_data(n[0:i], I[0:i])
        lsigma.set_data(n[0:i], sigma[0:i])
            
    
    anim = animation.FuncAnimation(fig, animate, init_func=init, fargs=(I,sigma))
    
    return anim



"""
And we can even use a weighted MC method, with weight function
   w = 1 - 0.64 x
that has to be normalized.
"""

#Method that performs the standard rectangle integration
def weighted_MC_method_integration(integrand, interval, N, samples, squaredsamples, weight, x):
    """
        integrand: integrand function
        interval: tuple containing the extremes (a,b) of the integration interval
        N: number of random numbers generated
        samples, squaredsamples: current available samples and squared of the samples
        
    Routine that performs the numerical integration of an integrand function in a certain interval (a,b) with the
    Monte Carlo method.
    Return the tuple (I, sigma, samplesn, squaredsamplesn), where:
        - I is the current value of the integral
        - sigma is the current integration error
        - samplesn and squaredsamplesn are the updated list of samples and their squares
    """
    
    #YOUR CODE GOES HERE





def execute_part5():


    #Define the parameters of the integration
    Npercycle = #YOUR CODE GOES HERE
    Ncycles = #YOUR CODE GOES HERE
    a = #YOUR CODE GOES HERE
    b = #YOUR CODE GOES HERE
    
    #Define the integrand function
    def integrand(x):
        #YOUR CODE GOES HERE
    
    #Define the weight function (normalized)
    def weight(x):
        #YOUR CODE GOES HERE
    
    #Define and the x(y) tranformation
    def x(y):
        #YOUR CODE GOES HERE
    
    
    n = np.array([]) #iterations
    I = np.array([]) #integrand
    sigma = np.array([]) #integration error
    
    samples = np.array([]) #array storing the samples
    squaredsamples = np.array([])
    
    for i in range(Ncycles):
        
        (In, sigman, samplesn, squaredsamplesn) = weighted_MC_method_integration(integrand, (a,b), Npercycle, samples, squaredsamples, weight, x)
    
        n = np.append(n,i)
        I = np.append(I,In)
        sigma = np.append(sigma,sigman)
    
        samples = samplesn
        squaredsamples = squaredsamplesn
 
    
    #First set up the figure, the axis, and the plot element we want to animate
    stdsize = plt.rcParams.get('figure.figsize')
    fig = plt.figure(figsize=(2*stdsize[0],stdsize[1]))
    axI = fig.add_subplot(121)
    axI.set_xlim(0, Ncycles)
    axI.set_ylim(0.60,0.65)
    lI, = axI.plot([],[], 'b-')
    axsigma = fig.add_subplot(122)
    axsigma.set_xlim(0, Ncycles)
    axsigma.set_ylim(0.,1.)
    lsigma, = axsigma.plot([],[],'r-')
    # YOU CAN ADJUST THE LABELS / LIMITS OF THE FIGURE YOURSELF
    
    
    def init():
        lI.set_data([], [])
        lsigma.set_data([],[])
    
    #And the animation function which is called every frame
    def animate(i, I, sigma):

        lI.set_data(n[0:i], I[0:i])
        lsigma.set_data(n[0:i], sigma[0:i])
            
    
    anim = animation.FuncAnimation(fig, animate, init_func=init, fargs=(I,sigma))    
    
    return anim

