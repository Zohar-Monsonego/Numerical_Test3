from numpy.linalg import norm

from gauss_seidel_formula import *


def find_seidel_g(A):
    _L = l_mat(A)  # calculate the L matrix
    _D = diagonal_mat(A)  # calculate the D matrix
    _U = U_mat(A)  # calculate the U matrix

    sum_L_D = sum_matrices(_L, _D)  # calculate: (L+D)
    inverse_sum_L_D = inverse(sum_L_D)  # calculate: (L+D)^-1
    minus_inverse_sum_L_D = mult_matrix_in_scalar(inverse_sum_L_D, -1)  # calculate: -(L+D)^-1
    minus_mult_U = matrix_multiply(minus_inverse_sum_L_D, _U)  # calculate: -(L+D)^-1 * U
    return minus_mult_U


def gauss_seidel(A, b, X0, TOL=1e-16, N=200):
    global x_rounded
    n = len(A)
    k = 1
    norm_gauss_g = MaxNorm(find_seidel_g(A))
    if norm_gauss_g > 1:
        print("G's norm: " + str(norm_gauss_g), "\n")
        print("The G's norm is too high, try doing a full pivot to reduce it")
    else:
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
        return tuple(x_rounded), norm_gauss_g


if __name__ == '__main__':
    A = np.array([[2, 3, 4, 5, 6], [-5, 3, 4, -2, 3], [4, -5, -2, 2, 6], [4, 5, -1, -2, -3], [5, 5, 3, -3, 5]])
    b = np.array([70, 20, 26, -12, 37])
    # X0 = np.zeros_like(b)
    B = np.array(([9, 3, 1], [4, 2, 1], [1, 1, 1]))
    b1 = np.array([-1, 4, 3])

    b_b = np.array([[9, 3, 1],
                    [4, 2, 1],
                    [1, 1, 1]])
    vector = np.array([-1, 4, 3])
    X0 = np.zeros_like(vector)

    A_ = np.array([[6, 3, 2, 4, 5], [-3, 5, 4, -1, -2], [5, 5, 5, 3, -3], [3, 3, -5, 4, -2], [6, -5, 4, -2, 2]])
    b_ = np.array([70, -12, 37, 20, 26])

    # X0 = np.zeros_like(b_)

    H = np.array([[5, 5, 5, -3, 3], [4, 5, -3, -2, -1], [4, -5, 6, 2, -2], [2, 3, 6, 5, 4], [-5, 3, 3, -2, 4]])
    h_ = np.array([37, -12, 26, 70, 20])
    # X0 = np.zeros_like(h_)

    F = np.array([[5, 5, -3, 5, 3], [4, 5, -2, -3, -1], [2, 3, 5, 6, 4], [4, -5, 2, 6, -2], [-5, 3, -2, 3, 4]])
    f_ = np.array([37, -12, 70, 26, 20])
    # X0 = np.zeros_like(f_)

    solution, g_norm = gauss_seidel(b_b, vector, X0)

    print(bcolors.OKBLUE, "\nApproximate solution:", solution)
    print("G's norm: " + str(g_norm))
