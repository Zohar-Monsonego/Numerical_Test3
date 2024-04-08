from colors import bcolors
from matrix_utility import row_addition_elementary_matrix, scalar_multiplication_elementary_matrix
import numpy as np
from matrix_utility import partial_pivoting


"""
Function that find the inverse of non-singular matrix
The function performs elementary row operations to transform it into the identity matrix. 
The resulting identity matrix will be the inverse of the input matrix if it is non-singular.
 If the input matrix is singular (i.e., its diagonal elements become zero during row operations), it raises an error.
"""

def inverse(matrix):
    print(bcolors.OKBLUE, f"=================== Finding the inverse of a non-singular matrix using elementary row operations ===================\n {matrix}\n", bcolors.ENDC)
    if matrix.shape[0] != matrix.shape[1]:  # checks if the row's number is equal to the col's number
        raise ValueError("Input matrix must be square.")
    else:
        n = matrix.shape[0]  # puts the row's number into a parameter
        identity = np.identity(n)  # creating an identity matrix of size n by n
        # Perform row operations to transform the input matrix into the identity matrix

        #for i in range(n):
        #if matrix[i, i] == 0:  # checks if the hinge number in the matrix equal 0
        #raise ValueError("Matrix is singular, cannot find its inverse.")

    for i in range(3):
        if matrix[i, i] == 0:
            print("We do pivoting ===")
            matrix, identity = partial_pivoting(matrix, i, n, identity)
            print(bcolors.OKGREEN,
            "------------------------------------------------------------------------------------------------------------------",
            bcolors.ENDC)
            print("-"
                  "-"
                  "-"
                  "identity matrix after pivoting: \n" + str(identity))

            print(bcolors.OKGREEN,
            "------------------------------------------------------------------------------------------------------------------",
            bcolors.ENDC)

        if matrix[i, i] != 1:
        # Scale the current row to make the diagonal element 1
            scalar = 1.0 / matrix[i, i]  # divide the hinge number by 1 to find the opposite number(example: 4-> -4)
            elementary_matrix = scalar_multiplication_elementary_matrix(n, i, scalar)
            print(f"elementary matrix to make the diagonal element 1 :\n {elementary_matrix} \n")
            matrix = np.dot(elementary_matrix, matrix)  # multiply the original matrix with the elementary matrix
            print(f"The matrix after elementary operation :\n {matrix}")
            print(bcolors.OKGREEN, "------------------------------------------------------------------------------------------------------------------",  bcolors.ENDC)
            identity = np.dot(elementary_matrix, identity)

        # Zero out the elements above and below the diagonal
        for j in range(n):
            if i != j:  # checks if the number is not a hinge number
                scalar = -matrix[j, i]   # make the minus of this number   ++++++=
                elementary_matrix = row_addition_elementary_matrix(n, j, i, scalar)
                print(f"elementary matrix for R{j+1} = R{j+1} + ({scalar}R{i+1}):\n {elementary_matrix} \n")
                matrix = np.dot(elementary_matrix, matrix)
                print(f"The matrix after elementary operation :\n {matrix}")
                print(bcolors.OKGREEN, "------------------------------------------------------------------------------------------------------------------",
                        bcolors.ENDC)
                identity = np.dot(elementary_matrix, identity)
                print("identity matrix: \n" + str(identity))
                print(bcolors.OKGREEN,
                        "------------------------------------------------------------------------------------------------------------------",
                        bcolors.ENDC)

    return identity


if __name__ == '__main__':

    A = np.array([[8, 1, 9],
                  [0, 0, 3],
                  [1, 0, 0]])

    B = np.array([[-1, 1, 3, -3, 1],
                  [3, -3, -4, 2, 3],
                  [2, 1, -5, -3, 5],
                  [-5, -6, 4, 1, 3],
                  [3, -2, -2, -3, 5]])

    try:
        A_inverse = inverse(B)
        print(bcolors.OKBLUE, "\nInverse of matrix A: \n", A_inverse)
        print("=====================================================================================================================", bcolors.ENDC)

    except ValueError as e:
        print(str(e))


