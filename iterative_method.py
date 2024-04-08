import math
import numpy as np
from colors import bcolors
import random
import math


# class bcolors:
#     OKBLUE = '\033[94m'
#     ENDC = '\033[0m'

def iterative_method(f, x0, lower_bound, upper_bound, tol=1e-6):
    """
    Implements the fixed-point iteration method to find roots of the equation f(x) = 0 within a specified range.

    Args:
    f (function): The function for which roots are to be found.
    x0 (float): Initial guess.
    lower_bound (float): Lower bound of the search range.
    upper_bound (float): Upper bound of the search range.
    tol (float): Tolerance (stop when |x_{n+1} - x_{n}| < tol).

    Returns:
    float: Approximation of the root.
    int: Number of iterations performed.
    """
    # print(f"Iterating with x0 = {x0} within range [{lower_bound}, {upper_bound}]")
    # print("Iter    xr        g(xr)      xr+1")
    xr = x0
    for n in range(10000):  # Iterating for 10 times
        g_xr = f(xr)
        xr_plus_1 = g_xr + xr

        # print(f"{n+1:<5}   {xr:<10.6f} {g_xr:<10.6f} {xr_plus_1:<10.6f}")  # Printing iteration
        if abs(xr_plus_1 - xr) < tol:
            return xr_plus_1, n

        xr = xr_plus_1
        if xr_plus_1 < lower_bound or xr_plus_1 > upper_bound:
            # print("Root is outside the specified range.")
            return None, -1

    # print("Series did not converge within the specified range.")
    return None, -1  # Return -1 to indicate non-convergence


def find_roots_iterative_method(a, b, c):
    f = lambda x: math.exp(x) - 1.5 * x - 1.6
    f1 = lambda x: a * x ** 2 + b * x ** 3 + c ** 2 * x ** 9
    f2 = lambda x: a * x ** 2 - b * math.sin(x) + c
    f3 = lambda x: math.sin(x)
    f4 = lambda x: x ** 2

    for i in range(-50, 50, 1):
        initial_guess = i
        root, iterations = iterative_method(f4, initial_guess, -50, 50)  # Adjust lower and upper bounds as needed
        if root is not None:
            print(bcolors.OKBLUE, f"Root:", root, bcolors.ENDC)
            # print(bcolors.OKBLUE,"Number of iterations:", iterations,bcolors.ENDC)
        # else:
        #     print(bcolors.OKBLUE, "Method did not converge within the specified range.", bcolors.ENDC)
        if i == 100 and root is None:
            print(bcolors.OKBLUE, "Method did not converge within the specified range.", bcolors.ENDC)


if __name__ == '__main__':
    find_roots_iterative_method(1, 0, 0)
