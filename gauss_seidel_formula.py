import numpy as np
from colors import bcolors
from matrix_utility import *


def inverse(matrix):
    # print(bcolors.OKBLUE, f"=================== Finding the inverse of a non-singular matrix using elementary row operations ===================\n {matrix}\n", bcolors.ENDC)
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Input matrix must be square.")

    n = matrix.shape[0]
    identity = np.identity(n)

    # Perform row operations to transform the input matrix into the identity matrix
    for i in range(n):
        if matrix[i, i] == 0:
            raise ValueError("Matrix is singular, cannot find its inverse.")

        if matrix[i, i] != 1:
            # Scale the current row to make the diagonal element 1
            scalar = 1.0 / matrix[i, i]
            elementary_matrix = scalar_multiplication_elementary_matrix(n, i, scalar)
            # print(f"elementary matrix to make the diagonal element 1 :\n {elementary_matrix} \n")
            matrix = np.dot(elementary_matrix, matrix)
            # print(f"The matrix after elementary operation :\n {matrix}")
            # print(bcolors.OKGREEN, "------------------------------------------------------------------------------------------------------------------",  bcolors.ENDC)
            identity = np.dot(elementary_matrix, identity)

        # Zero out the elements above and below the diagonal
        for j in range(n):
            if i != j:
                scalar = -matrix[j, i]
                elementary_matrix = row_addition_elementary_matrix(n, j, i, scalar)
                # print(f"elementary matrix for R{j+1} = R{j+1} + ({scalar}R{i+1}):\n {elementary_matrix} \n")
                matrix = np.dot(elementary_matrix, matrix)
                # print(f"The matrix after elementary operation :\n {matrix}")
                # print(bcolors.OKGREEN, "------------------------------------------------------------------------------------------------------------------",
                #       bcolors.ENDC)
                identity = np.dot(elementary_matrix, identity)

    return identity


def UMatrix(matrix, vector):
    """
    :param matrix: Matrix nxn
    :return:Disassembly into a  U matrix
    """
    # result matrix initialized as singularity matrix
    U = MakeIMatrix(len(matrix), len(matrix))
    # loop for each row
    for i in range(len(matrix[0])):
        # pivoting process
        matrix, vector = RowXchageZero(matrix, vector)
        for j in range(i + 1, len(matrix)):
            elementary = MakeIMatrix(len(matrix[0]), len(matrix))
            # Finding the M(ij) to reset the organs under the pivot
            elementary[j][i] = -(matrix[j][i]) / matrix[i][i]
            matrix = MultiplyMatrix(elementary, matrix)
    # U matrix is a doubling of elementary matrices that we used to reset organs under the pivot
    U = MultiplyMatrix(U, matrix)
    return U


def LMatrix(matrix, vector):
    """
       :param matrix: Matrix nxn
       :return:Disassembly into a  L matrix
       """
    # Initialize the result matrix
    L = MakeIMatrix(len(matrix), len(matrix))
    # loop for each row
    for i in range(len(matrix[0])):
        # pivoting process
        matrix, vector = RowXchageZero(matrix, vector)
        for j in range(i + 1, len(matrix)):
            elementary = MakeIMatrix(len(matrix[0]), len(matrix))
            # Finding the M(ij) to reset the organs under the pivot
            elementary[j][i] = -(matrix[j][i]) / matrix[i][i]
            # L matrix is a doubling of inverse elementary matrices
            L[j][i] = (matrix[j][i]) / matrix[i][i]
            matrix = MultiplyMatrix(elementary, matrix)

    return L


def l_mat(A):
    rows, cols = A.shape
    B = np.zeros_like(A)

    for i in range(rows):
        for j in range(cols):
            if i > j:  # This condition checks if the element is below the diagonal
                B[i, j] = A[i, j]

    return B


def U_mat(A):
    rows, cols = A.shape
    B = np.zeros_like(A)

    for i in range(rows):
        for j in range(cols):
            if i < j:  # This condition checks if the element is below the diagonal
                B[i, j] = A[i, j]

    return B


def diagonal_mat(A):
    rows, cols = A.shape
    B = np.zeros_like(A)

    for i in range(rows):
        for j in range(cols):
            if i == j:  # This condition checks if the element is below the diagonal
                B[i, j] = A[i, j]

    return B


def sum_matrices(A, B):
    size_A = A.shape
    size_B = B.shape
    rows, cols = size_A
    if not size_A == size_B:
        print("can not perform operation, sizes are not equal")
        return

    sum_mat = np.zeros_like(A)

    for i in range(rows):
        for j in range(cols):
            sum_mat[i][j] = A[i][j] + B[i][j]
    return sum_mat


def sum_vectors(a, b):
    if isinstance(a, list):
        size_A = len(a)
        size_B = len(b)
    elif isinstance(a, np.ndarray):
        size_A = a.shape[0]
        size_B = b.shape[0]
    else:
        raise TypeError("Unsupported type for 'a'. Must be a list or a NumPy array.")
    if not size_A == size_B:
        print("can not perform operation, sizes are not equal")
        return
    rows = size_A

    sum_v = np.zeros_like(a)

    for i in range(rows):
        sum_v[i][0] = a[i][0] + b[i][0]
    return sum_v


def mult_matrix_in_scalar(A, scalar):
    size_A = A.shape
    rows, cols = size_A

    mult_mat = np.zeros_like(A)

    for i in range(rows):
        for j in range(cols):
            mult_mat[i][j] = A[i][j] * scalar
    return mult_mat



def MulMatrixVector(InversedMat, b_vector):
    """
    Function for multiplying a vector matrix
    :param InversedMat: Matrix nxn
    :param b_vector: Vector n
    :return: Result vector
    """
    result = []
    # Initialize the x vector
    for i in range(len(b_vector)):
        result.append([])
        result[i].append(0)
    # Multiplication of inverse matrix in the result vector
    for i in range(len(InversedMat)):
        for k in range(len(b_vector)):
            result[i][0] += InversedMat[i][k] * b_vector[k][0]
    return result


def gauss_seidel_formula(matrix, Xr, b):
    _L = l_mat(matrix)  # calculate the L matrix
    _D = diagonal_mat(matrix)  # calculate the D matrix
    _U = U_mat(matrix)  # calculate the U matrix

    sum_L_D = sum_matrices(_L, _D)  # calculate: (L+D)
    inverse_sum_L_D = inverse(sum_L_D)  # calculate: (L+D)^-1
    inverse_mult_U = matrix_multiply(inverse_sum_L_D, _U)  # calculate: (L+D)^-1 * U
    becomes_minus = mult_matrix_in_scalar(inverse_mult_U, -1)  # calculate: -(L+D)^-1 * U
    mat_mult_x0 = MulMatrixVector(becomes_minus, Xr)  # calculate: (-(L+D)^-1 * U) *X0

    inverse_mult_b = MulMatrixVector(inverse_sum_L_D, b)  # calculate: (L+D)^-1 * b

    result_vector = sum_vectors(mat_mult_x0, inverse_mult_b)
    return result_vector


if __name__ == '__main__':
    A = np.array([[5, 1, 2],
                  [1, 6, 4],
                  [0, 3, 8]])

    b = np.array([[1],
                  [2],
                  [3]])

    Xr = np.array([[0],
                   [0],
                   [0]])

    if is_diagonally_dominant(A):
        while True:
            Xr_plus_1 = np.array(gauss_seidel_formula(A, Xr, b))
            if not np.any(np.abs(Xr - Xr_plus_1) > 0.001):
                print(Xr_plus_1)
                break
            Xr = Xr_plus_1.copy()
    else:
        print("The diagonal in the matrix is not dominant and therefore the iterative methods will not converge")
