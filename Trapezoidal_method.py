import math

from colors import bcolors


def trapezoidal_rule(f, a, b, n):
    h = (b - a) / n
    T = f(a) + f(b)
    integral = 0.5 * T  # Initialize with endpoints

    for i in range(1, n):
        x_i = a + i * h
        integral += f(x_i)

    integral *= h

    return integral


def calculate_error(e, n, f):
    global new_result, correct_result
    i = n
    while i < 200:
        result1 = trapezoidal_rule(f, 0, 1, i)
        result2 = trapezoidal_rule(f, 0, 1, i + 1)
        print("result 1 = " + str(result1))
        print("result 2 = " + str(result2))
        while abs(result2 - result1) > e:
            difference = abs(result2 - result1)
            print("difference = " + str(difference))
            result1 = result2
            i += 1
            result2 = trapezoidal_rule(f, 0, 1, i + 1)
            print(bcolors.OKGREEN, "New result = " + str(result2), bcolors.ENDC)
        correct_result = result2
        break

    return correct_result


if __name__ == '__main__':
    f = lambda x: math.e ** (x ** 2)
    error = 0.002
    n = 2
    result = calculate_error(error, n, f)

    print(bcolors.OKBLUE, "Approximate integral:", result, bcolors.ENDC)

