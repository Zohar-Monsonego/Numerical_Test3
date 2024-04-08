import numpy as np
from numpy.linalg import norm

from gauss_seidel_formula import *
from Matrix.inverse_matrix import inverse

from colors import bcolors


def find_jacobi_g(A):
    L = l_mat(A)  # calculate the L matrix
    D = diagonal_mat(A)  # calculate the D matrix
    U = U_mat(A)  # calculate the U matrix
    L_plus_U = sum_matrices(L, U)  # calculate: (L+U)
    invesre_L_plus_U = inverse(L_plus_U)  # calculate : (L+U)^-1
    inverse_D = inverse(D)  # calculate : D^-1
    minus_inverse_D = mult_matrix_in_scalar(inverse_D, -1)  # calculate: -D^-1
    minus_inverse_D_mult_inverse_L_plus_U = matrix_multiply(minus_inverse_D,
                                                            invesre_L_plus_U)  # calculate : -D^-1 * (L+U)^-1
    return minus_inverse_D_mult_inverse_L_plus_U


def jacobi_iterative(A, b, X0, TOL=1e-16, N=250):
    global x_rounded
    n = len(A)
    k = 1
    norm_jacobi_g = norm(find_jacobi_g(A))
    # if norm_jacobi_g > 1:
    # print("G's norm: " + str(norm_jacobi_g), "\n")
    # print("The G's norm is too high, try doing a full pivot to reduce it")
    # else:
    print(
        "Iteration" + "\t\t\t".join(
            [" {:>12}".format(var) for var in ["x{}".format(i) for i in range(1, len(A) + 1)]]))
    print("-----------------------------------------------------------------------------------------------")
    while k <= N:
        x = np.zeros(n, dtype=np.double)
        for i in range(n):
            sigma = 0
            for j in range(n):
                if j != i:
                    sigma += A[i][j] * X0[j]
            x[i] = (b[i] - sigma) / A[i][i]
        x_rounded = np.round(x)

        print("{:<15} ".format(k) + "\t\t".join(["{:<15} ".format(val) for val in x]))

        if norm(x - X0, np.inf) < TOL:
            return tuple(x_rounded)

        k += 1
        X0 = x.copy()

    print("Maximum number of iterations exceeded")
    return tuple(x_rounded)


if __name__ == "__main__":
    A = np.array([[3, -1, 1], [0, 1, -1], [1, 1, -2]])  # here, G = 1.6 and it does converge
    b = np.array([4, -1, -3])

    x1 = np.zeros_like(b, dtype=np.double)

    b_b = np.array([[9, 3, 1],  # in this matrix, the G = 1.8 and it doesn't converge
                    [4, 2, 1],
                    [1, 1, 1]])
    vector = np.array([-1, 4, 3])
    X0 = np.zeros_like(vector, dtype=np.double)

    solution = jacobi_iterative(A, b, x1)

    norm = MaxNorm(find_jacobi_g(A))
    print("norm: " + str(norm))

    print(bcolors.OKBLUE, "\nApproximate solution:", solution)
