from colors import bcolors
from matrix_utility import *
from gauss_seidel2 import gauss_seidel


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


def SolveLU(matrix, vector):
    """
    Function for deconstructing a linear equation by ungrouping LU
    :param matrix: Matrix nxn
    :param vector: Vector n
    :return: Solve Ax=b -> x=U(-1)L(-1)b
    """
    matrixU = UMatrix(matrix)
    matrixL = LMatrix(matrix)
    return MultiplyMatrix(InverseMatrix(matrixU), MultiplyMatrix(InverseMatrix(matrixL), vector))


def solveMatrix(matrixA, vectorb):  # changed the solution to gauss seidel method and not gauss jordan elimination
    detA = Determinant(matrixA, 1)
    print(bcolors.YELLOW, "\nDET(A) = ", detA)

    if detA != 0:
        print(bcolors.OKBLUE, "\nnon-Singular Matrix - Perform Gauss Seidel method", bcolors.ENDC)
        result = solve_gauss_seidel(matrixA, vectorb)
        print(np.array(result))
        return result
    else:
        print("Singular Matrix - Perform LU Decomposition\n")
        print("Matrix U: \n")
        print(np.array(UMatrix(matrixA, vectorb)))
        print("\nMatrix L: \n")
        print(np.array(LMatrix(matrixA, vectorb)))
        print("\nMatrix A=LU: \n")
        result = MultiplyMatrix(LMatrix(matrixA, vectorb), UMatrix(matrixA, vectorb))
        print(np.array(result))
        return result


def polynomialInterpolation(table_points, x):
    matrix = [[point[0] ** i for i in range(len(table_points))] for point in table_points]  # Makes the initial matrix

    b = [[point[1]] for point in table_points]

    print(bcolors.OKBLUE, "The matrix obtained from the points: ", bcolors.ENDC, '\n', np.array(matrix))
    print(bcolors.OKBLUE, "\nb vector: ", bcolors.ENDC, '\n', np.array(b))
    matrixSol = solveMatrix(matrix, b)

    result = sum([matrixSol[i] * (x ** i) for i in range(len(matrixSol))])
    print(bcolors.OKBLUE, "\nThe polynom:", bcolors.ENDC)
    print('P(X) = ' + '+'.join(['(' + str(matrixSol[i]) + ') * x^' + str(i) + ' ' for i in range(len(matrixSol))]))
    print(bcolors.OKGREEN, f"\nThe Result of P(X={x}) is:", bcolors.ENDC)
    print(result)
    return result


def solve_gauss_seidel(A, b):
    x0 = np.zeros_like(b)
    result = gauss_seidel(A, b, x0)
    return result


if __name__ == '__main__':
    table_points = [(0, 0), (1, 0.8415), (2, 0.9093), (3, 0.1411), (4, -0.7568), (5, -0.9589), (6, -0.2794)]
    x = 1.28
    # table_points = [(1, 3), (2, 4), (3, -1)]
    # x = 1.5
    print(bcolors.OKBLUE, "----------------- Interpolation & Extrapolation Methods -----------------\n", bcolors.ENDC)
    print(bcolors.OKBLUE, "Table Points: ", bcolors.ENDC, table_points)
    print(bcolors.OKBLUE, "Finding an approximation to the point: ", bcolors.ENDC, x, '\n')
    polynomialInterpolation(table_points, x)
    print(bcolors.OKBLUE, "\n---------------------------------------------------------------------------\n",
          bcolors.ENDC)
