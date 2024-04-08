from colors import bcolors


def newton_raphson(f, df, p0, TOL, N=50):
    print("{:<10} {:<15} {:<15} ".format("Iteration", "po", "p1"))
    for i in range(N):
        if df(p0) == 0:
            print("Derivative is zero at p0, method cannot continue.")
            return

        p = p0 - f(p0) / df(p0)

        if abs(p - p0) < TOL:
            return p  # Procedure completed successfully
        print("{:<10} {:<15.9f} {:<15.9f} ".format(i, p0, p))
        p0 = p
    return p


if __name__ == '__main__':
    f = lambda x: x**2 - 5*x + 2
    df = lambda x: 2*x - 5
    p0 = 4  # first guess, everytime we need to change the first guess and to check if there is another roots
    TOL = 1e-6
    N = 100
    roots = newton_raphson(f, df,p0,TOL,N)
    print(bcolors.OKBLUE,"\nThe equation f(x) has an approximate root at x = {:<15.9f} ".format(roots),bcolors.ENDC,)
