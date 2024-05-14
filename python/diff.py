import numpy as np

def step(f, t, x, h): 
    kn1 = f(t, x)
    kn2 = f(t + h/2, x + (h/2)*kn1)
    kn3 = f(t + h/2, x + (h/2)*kn2)
    kn4 = f(t + h, x + h*kn3)

    return x + (h/6)*(kn1 + 2*kn2 + 2*kn3 + kn4)

def solve_diffeq(f, t0, x0, h, maxiter, terminate):
    n = len(x0)
    tn = t0
    xn = x0

    sol = np.array([t0, *x0]).reshape(1, n+1)

    for i in range(1, maxiter):
        next_x = step(f, tn, xn, h)
        tn += h
        if terminate(i, tn, xn):
            break

        xn = next_x
        sol = np.append(sol, np.array([tn, *xn]).reshape(1, n+1), axis=0)

    return sol