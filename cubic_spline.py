import jacobi_utilities
from sympy import *

x = Symbol('x')


def natural_cubic_spline(f, x0):
    global s
    h = list()
    for i in range(len(f) - 1):
        h.append(f[i + 1][0] - f[i][0])  # save the h for each two points

    m = list()  # calculating: (1/6) * h[i - 1]
    m.append(0)
    for i in range(1, len(f) - 1):
        m.append((1 / 6) * h[i - 1])

    z = list()  # calculating: (1/3) * (h[i-1] + h[i])
    z.append(0)
    for i in range(1, len(f) - 1):
        z.append((1 / 3) * (h[i - 1] + h[i]))

    w = list()  # calculating: (1/6) * h[i]
    w.append(0)
    for i in range(1, len(f) - 1):
        w.append((1 / 6) * h[i])

    d = list()
    d.append(0)  # d0=0
    for i in range(1, len(f) - 1):
        d.append((f[i + 1][1] - f[i][1]) / h[i] - (f[i][1] - f[i - 1][1]) / h[i - 1])  # calculating the vector b
    d.append(0)  # dn

    # building the matrix
    mat = list()

    # first row
    mat.append(list())
    mat[0].append(1)
    for j in range(len(f) - 1):
        mat[0].append(0)

    for i in range(1, len(f) - 1):
        mat.append(list())
        for j in range(len(f)):
            if j == i - 1:  # put miu
                mat[i].append(m[i])
            elif j == i:
                mat[i].append(z[i])
            elif j == i + 1:  # put lambda
                mat[i].append(w[i])
            else:
                mat[i].append(0)

    # last row
    mat.append(list())
    for j in range(len(f) - 1):
        mat[len(f) - 1].append(0)
    mat[len(f) - 1].append(1)

    print("matrix: " + str(mat))
    print("vector b: " + str(d))

    # get m vector
    print("\nJacobi middle results: ")
    M = (jacobi_utilities.Jacobi(mat, d))
    print("\nvector M: " + str(list(map(float, M))))

    # find S:
    for i in range(1, len(f)):
        s = (((f[i][1]) * (x - f[i - 1][0])) / h[i - 1]) - ((f[i - 1][1]) * (x - f[i][0])) / (h[i - 1])
        s += (((M[i - 1]) / 6) * ((((x - f[i][0]) ** 3) / (h[i - 1])) - ((h[i - 1]) * (x - f[i][0]))))
        s -= ((M[i]) / 6) * ((((x - f[i - 1][0]) ** 3) / (h[i - 1])) - h[i - 1] * (x - f[i - 1][0]))
        print("s" + str(i - 1) + "(x) = " + str(s))

    # find the location of x0:
    loc = 0
    for i in range(1, len(f)):
        if x0 < f[i][0] and x0 > f[i - 1][0]:
            loc = i
            print(
                "\nx0 between f(x" + str(loc - 1) + ") = " + str(f[loc - 1][0]) + " and f(x" + str(loc) + ") = " + str(
                    f[loc][0]) + " so:")
            print("s" + str(loc - 1) + "(" + str(x0) + ") = " + str(float(s.subs(x, x0))))
            break

        if loc == 0:
            print("no range found for x0")
            return


if __name__ == '__main__':
    f = [(1, 1), (2, 2), (3, 1), (4, 1.5), (5, 1)]
    x0 = 1.5

    print("func: " + str(f))
    print("x0 = " + str(x0) + "\n")
    natural_cubic_spline(f, x0)

