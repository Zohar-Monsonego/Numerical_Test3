from numpy.linalg import norm
import numpy as np
from colors import bcolors
#from gauss_seidel_formula import *


def gauss_seidel(A, b, X0, TOL=1e-16, N=200):
    global x_rounded
    n = len(A)
    k = 1

    print(
        "Iteration" + "\t\t\t".join(
            [" {:>12}".format(var) for var in ["x{}".format(i) for i in range(1, len(A) + 1)]]))
    print("-----------------------------------------------------------------------------------------------")
    x = np.zeros(n, dtype=np.double)
    while k <= N:

        for i in range(n):
            sigma = 0
            for j in range(n):
                if j != i:
                    sigma += A[i][j] * x[j]
            x[i] = (b[i] - sigma) / A[i][i]
        x_rounded = np.round(x)

        print("{:<15} ".format(k) + "\t\t".join(["{:<15} ".format(val) for val in x]))

        if norm(x - X0, np.inf) < TOL:
            return tuple(x_rounded)

        k += 1
        X0 = x.copy()

    print("Maximum number of iterations exceeded")
    return tuple(x_rounded)


if __name__ == '__main__':
    A = np.array([[2, 3, 4, 5, 6], [-5, 3, 4, -2, 3], [4, -5, -2, 2, 6], [4, 5, -1, -2, -3], [5, 5, 3, -3, 5]])
    b = np.array([70, 20, 26, -12, 37])
    #X0 = np.zeros_like(b)
    B = np.array(([9, 3, 1], [4, 2, 1], [1, 1, 1]))
    b1 = np.array([-1, 4, 3])

    b_b = np.array([[9, 3, 1],
                    [4, 2, 1],
                    [1, 1, 1]])
    vector = np.array([-1, 4, 3])
    X0 = np.zeros_like(vector)

    solution = gauss_seidel(b_b, vector, X0)

    print(bcolors.OKBLUE, "\nApproximate solution:", solution)

