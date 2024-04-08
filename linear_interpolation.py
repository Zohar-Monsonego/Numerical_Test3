from colors import bcolors


def linearInterpolation(table_points, point):
    for i in range(len(table_points) - 1):
        if i <= point <= i + 1:
            x1 = table_points[i][0]  # insert into x1 the y of the i point
            x2 = table_points[i + 1][0]
            y1 = table_points[i][1]
            y2 = table_points[i + 1][1]
            result = (((y1 - y2) / (x1 - x2)) * point) + ((y2 * x1) - (y1 * x2)) / (x1 - x2)
            print(bcolors.OKGREEN, "\nThe approximation (interpolation) of the point ", point, " is: ",
                  bcolors.ENDC,
                  round(result, 4))


if __name__ == '__main__':
    #table_points = [(0, 0), (1, 0.8415), (2, 0.9093), (3, 0.1411), (4, -0.7568), (5, -0.9589), (6, -0.2794)]
    #x = 1.28
    table_points = [(1, 3), (2, 4), (3, -1)]
    x = 1.5
    print(bcolors.OKBLUE, "----------------- Interpolation & Extrapolation Methods -----------------\n", bcolors.ENDC)
    print(bcolors.OKBLUE, "Table Points: ", bcolors.ENDC, table_points)
    print(bcolors.OKBLUE, "Finding an approximation to the point: ", bcolors.ENDC, x)
    linearInterpolation(table_points, x)
    print(bcolors.OKBLUE, "\n---------------------------------------------------------------------------\n",
          bcolors.ENDC)

