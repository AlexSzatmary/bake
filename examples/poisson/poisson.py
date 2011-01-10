def thomas_solver(b):
    print b
    for i in xrange(1, len(b)):
        b[i] = b[i] + b[i - 1] * i/(i + 1.)
    print b
    x = [0. for i in xrange(len(b))]
    x[-1] = - b[-1] * (i + 1.)/(i + 2.)
    print range(len(b) - 1, -1, -1)
    for i in xrange(len(b) - 2, -1, -1):
        print i
        x[i] = -(b[i] - x[i + 1])*(i + 1)/(i + 2)
    return x

def set_up():
    x_l = 0.
    x_r = 1.
    u_l = 0.
    u_r = 0.
    n = 11
    h = (x_r - x_l)/(n - 1)
    def f(x):
        return 1.

    u_ideal = [(x_l + i*h)**2/2 - (x_l + i*h)/2 for i in xrange(0, n)]
    b = [f(x_l + h*i)*h*h for i in xrange(1, n - 1)]
    b[0] -= u_l
    b[-1] -= u_r
    u = thomas_solver(b)
    print(u_l)
    for u_i in u:
        print(u_i)
    print(u_r)

    print(u_l - u_ideal[0])
    for i in xrange(1, n - 1):
        print(u[i - 1] - u_ideal[i])
    print(u_r - u_ideal[-1])

    print (u_l - 2 * u[0] + u[1])/(h*h)
    for i in xrange(1, n - 3):
        print((u[i - 1] - 2*u[i] + u[i + 1])/(h*h))
    print (u[n - 4] - 2*u[n - 3] + u_r)/(h*h)
