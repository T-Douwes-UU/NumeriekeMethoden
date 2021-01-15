############# PART 1 ##############
"""
Let's consider a bullet (an angry bird) which is launched from a starting point (x0, y0) = (0,y0) with initial 
velocity (vx0,vy0) = (v0cos(theta),v0sin(theta)), namely with speed v0 and an angle theta in [0,90] degrees with 
respect to the ground.

Our aim is to determine the value of the launching angle theta to hit a target (a pig) placed in (xf,yf) = (xf, 0), 
for a given value of y0, v0 and xf.
Namely, if we call x(theta) the distance that the bird travels before landing for a certaing launching angle theta
(given the starting height y0 and the starting speed v0), we want to find the solution to the equation x(theta) - xf = 0.
From Newton's equations of motion we can determine
    x(theta) = (1./g) * [v0**2 * cos(theta) * sin(theta) + sqrt(v0**2 * sin(theta)**2 + 2 * y0 * g)]
where g is the gravity acceleration constant, namely g = 9.8 m / s**2.

If y0 is different from 0, the equation x(theta) - xf = 0 cannot be solved analytically, and we need to apply one of
the numerical methods discussed in the lectures.
However, in order to understand if a solution even exists for given values of y0, v0 and xf, and have an idea of its
position, we can plot the value of x(theta) - xf for the given parameters.
"""
#%%
import matplotlib.pyplot as plt
import numpy as np


# Define the value of g
g = 9.8  # m/s**2


# Define x(theta)
def xTheta(theta_deg, y0, v0):
    """
    theta_deg: launching angle in degrees
    y0: starting height
    v0: starting speed

    Returns the value of x(theta).
    """
    return v0*np.cos(theta_deg)/g * (v0*np.sin(theta_deg) + np.sqrt(v0**2 * np.sin(theta_deg)**2 + 2 * y0 * g))


def execute_part_1():
    # Get a set of 200 angles equally spaced within the interval [0,90] degrees
    angles = np.linspace(0, np.pi/2, 200)
    ticks = np.linspace(0, 1.5, 16)

    # Plot x(theta)-xf, with a line in zero as a reference for a set of parameters
    i=1
    for (y0, v0, xf) in [(0., 8., 10.2), (0., 10., 10.2), (0., 12., 10.2), (5., 5., 8), (5., 5., 5.664),
                         (5., 10., 12.), (8., 10., 12.)]:
        plt.figure(i)
        plt.title(r'y$_0$={0}, v$_0$ = {1}, x$_f$ = {2}'.format(y0, v0, xf))
        plt.plot(angles, np.array([(xTheta(theta, y0, v0)-xf) for theta in angles]))
        plt.plot(angles, np.zeros(200), 'r--')
        plt.xticks(ticks)
        plt.xlabel(r'$\theta$ (rad)')
        plt.ylabel(r'x($\theta$) - x$_f$')
        plt.savefig("plot{}.png".format(i), dpi=300, bbox_inches='tight')
        plt.show()
        i += 1


execute_part_1()

#%%

############# PART 2 ##############
"""
As can be seen from the plots, for different values of the parameters the equation x(theta)-xf = 0 has either 0, 1 or 2
roots.
This is true both when y0=0, when x(theta) is symmetrical around 45 degrees, or when y0 != 0 and x(theta) has a richer
and less trivial behaviour.
The plots allow us to determine intervals of thetas [thetai,thetai+1] which bracket the (eventual) solution(s). However,
the effective solution has to be determined numerically.
In order to determine the solution numerically we can for example use the bisection method.
"""

#Build the method
def bisection_method(function, interval, tolerance):
    """
    function: f(x) of which we want to find the root
    interval: a tuple (a,b) defining the starting interval which brackets the solution

    Function that determines the root of the equation f(x) = 0 within the interval (a, b) using the method
    of bisection.

    Returns a tuple (Root, Error, Number of Iterations)
    """
    a, b = interval
    error = b-a
    i = 0
    while error > tolerance:
        i += 1
        c = (a+b)/2.
        error *= 1/2
        if function(a)*function(c) < 0:
            b = c
        elif function(a)*function(c) > 0:
            a = c
        else:
            error = 0
            break
    return c, error, i


def execute_part_2():
    # Consider a set of parameters for which only one solution exists (given previous plots)
    y0, v0, xf = 8., 10., 12.

    # From the corresponding plot, see that the solution is located between
    theta0 = 0.9
    # and
    theta1 = 1.3

    # Try to find a solution in the interval using the bisection method
    def function(theta):
        return xTheta(theta, y0, v0) - xf
    interval = (theta0, theta1)
    print(bisection_method(function, interval, 1e-6))


execute_part_2()

#%%

############# PART 3 ##############
"""
The exact same problem can be faced using other methods, like the secant method.
Notice that in this case the starting interval [theta0,theta1] does not have to include the solution,
but it just has to be close to it.
"""

#Build the method
def secant_method(function, interval):
    """
    function: f(x) of which we want to find the root
    interval: a tuple (a,b) defining the starting interval, which in this case does not have to include the solution
              but just to be close to it

    Function that determines the root of the equation f(x) = 0 within the interval (a, b) using the secant method.

    Returns a tuple (Root, Error, Number of Iterations)
    """
    #YOUR CODE GOES HERE



def execute_part_3():
    #Consider a set of parameters for which only one solution exists (given previous plots)
    #(y0,v0,xf) = #YOUR CODE GOES HERE

    #From the corresponding plot, see that two points close to the solution are
    #theta0 = #YOUR CODE GOES HERE
    #and
    #theta1 = #YOUR CODE GOES HERE

    #Try to find a solution in the interval using the secant method
    #YOUR CODE GOES HERE



############# PART 4 ##############
    """
    You can also try to use different methods, as the Raphson-Newton method or the Brent's method! :)
    """

