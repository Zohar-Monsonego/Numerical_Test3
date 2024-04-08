import math
import numpy as np
import matplotlib.pyplot as plt

import sympy as sp

from colors import bcolors
from sympy.utilities.lambdify import lambdify

x = sp.symbols('x')


def simpsons_rule(f, a, b, n):
    if n % 2 != 0:
        raise ValueError("Number of subintervals (n) must be even for Simpson's Rule.")

    h = (b - a) / n

    integral = f(a) + f(b)  # Initialize with endpoints

    for i in range(1, n):
        x_i = a + i * h
        if i % 2 == 0:
            integral += 2 * f(x_i)
        else:
            integral += 4 * f(x_i)

    integral *= h / 3

    return integral


def calculate_error(e, n, f):
    global new_result, correct_result
    i = n
    while i < 200:
        result1 = simpsons_rule(f, 0, 1, i)
        result2 = simpsons_rule(f, 0, 1, i + 2)
        print("result 1 = " + str(result1))
        print("result 2 = " + str(result2))
        while abs(result2 - result1) > e:
            difference = abs(result2 - result1)
            print("difference = " + str(difference))
            result1 = result2
            i += 2
            result2 = simpsons_rule(f, 0, 1, i + 2)
            print(bcolors.OKGREEN, "New result = " + str(result2), bcolors.ENDC)
        correct_result = result2
        break

    return correct_result


if __name__ == '__main__':
    f = lambda x: math.e ** (x ** 2)
    n = 6
    a = 0
    b = 1
    error = 0.00001
    integral = calculate_error(error, n, f)

    print(f" Division into n={n} sections ")
    #integral = simpsons_rule(f, 0, 1, n)
    print(bcolors.OKBLUE, f"Numerical Integration of definite integral in range [{a},{b}] is {integral}", bcolors.ENDC)
