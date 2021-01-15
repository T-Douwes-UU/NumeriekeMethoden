import numpy as np
from numpy.random import default_rng
import matplotlib.pyplot as plt


rng = default_rng(2020)  # Seeded random number generator

n = 100  # Amount of points in the dataset
a = 1  # Slope of the linear dependence
std = 10  # Standard deviation of noise

x = np.arange(n) + 1  # Standard domain running from 1 to n


def generate_linear_data(domain=x, slope=a, noise_scale=std):
    """Generates a linear dataset with added random noise.

    Data is generated using the formula `y=a*x`, then an equal amount of normally distributed noise is added
    to each point.

    Args:
        domain (optional): The points x_i over which the data is generated. Defaults to x.
        slope (optional): The slope along which data is generated. Defaults to a.
        noise_scale (optional): The standard deviation of noise added to the clean data. Defaults to std.

    Returns:
        A 1-dimensional NumPy array containing the y data.

    """
    clean_data = slope * domain
    noise = rng.normal(scale=noise_scale, size=len(domain))
    data = clean_data + noise
    return data


y = generate_linear_data()  # The linear data set we'll be working with


def execute_part1():
    """Makes a scatter plot of our data set."""
    plt.plot(x, y, 'bo')
    plt.show()


def cost_function_gradient(w, b, data=y, domain=x):
    """Calculates the gradient of the cost function `J(w,b)` for a linear fitting function.

    Args:
        w: The slope of the linear fitting function.
        b: The constant term in the linear fitting function.
        data (optional): The data set that the function is fit to. Defaults to y.
        domain (optional): The domain of x_i over which the data resides. Defaults to x.

    Returns:
        A tuple consisting of the partial derivatives w.r.t. w and b, respectively.

    """
    predictions = w * domain + b
    diffs = predictions - data
    weighted_diffs = diffs * domain
    grad_w = sum(weighted_diffs)/len(data)
    grad_b = sum(diffs)/len(data)
    return grad_w, grad_b


def gradient_descent(w_start, b_start, rate, w_tolerance=1.e-6, b_tolerance=1.e-6,
                     data=y, domain=x):
    """Uses gradient descent to find a linear fit to a given dataset.

    Args:
        w_start: First guess for the slope w.
        b_start: First guess for the constant term b.
        rate: The learning rate for linear regression. Too high causes overshooting, too low makes the process slow.
            Typical values are quite low, around the order of 10^-4.
        w_tolerance (optional): The minimum desired precision for w. Defaults to 10^-6.
        b_tolerance (optional): The minimum desired precision for b. Defaults to 10^-6.
        data (optional): The data set that the function is fit to. Defaults to y.
        domain (optional): The domain of x_i over which the data resides. Defaults to x.

    Returns:
        A tuple containing w and b.
        Also prints a statement in the console containing w, b, and their tolerances.

    """
    w = w_start
    b = b_start
    w_precision = w_tolerance + 1
    b_precision = b_tolerance + 1
    count = 0
    while w_precision > w_tolerance or b_precision > b_tolerance:
        grad_w, grad_b = cost_function_gradient(w, b, data=data, domain=domain)
        w -= rate * grad_w
        b -= rate * grad_b
        w_precision = abs(grad_w)  # Not including the learning rate in the precision seems counterintuitive,
        b_precision = abs(grad_b)  # but this prevents early halting when the learning rate is low.
        count += 1
    print(f"Found w = {w}+/-{w_tolerance} and b = {b}+/-{b_tolerance} in {count} iterations")
    return w, b


def execute_part2():
    """Finds a best fit to our data set using linear regression, comparing it to a builtin NumPy result.
    The best fit is then plotted against the data set.
    """
    w, b = gradient_descent(1., 0., .0005)
    [w_np, b_np] = np.polyfit(x, y, 1)
    print(f"NumPy builtin np.polyfit(x,y,1) returns w = {w_np} and b = {b_np}")
    plt.plot(x, y, 'bo', x, w * x + b, 'r-')
    plt.show()
