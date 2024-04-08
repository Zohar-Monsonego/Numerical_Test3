from colors import bcolors


def lagrange_interpolation(x_data, y_data, x):
    n = len(x_data)
    result = 0.0

    for i in range(n):
        term = y_data[i]
        for j in range(n):
            if i != j:
                term *= (x - x_data[j]) / (x_data[i] - x_data[j])
        result += term

    return result


if __name__ == '__main__':
    x_data = [1, 2, 5]
    y_data = [1, 0, 2]

    x_interpolate = 3  # The x-value where you want to interpolate
    y_interpolate = lagrange_interpolation(x_data, y_data, x_interpolate)

    #x = [1, 2, 3]
    #y = [3, 4, -1]
    #x_ = 1.5
    #y_ = lagrange_interpolation(x, y, x_)
    print(bcolors.OKBLUE, "\nInterpolated value at x =", x_interpolate, "is y =", y_interpolate, bcolors.ENDC)
