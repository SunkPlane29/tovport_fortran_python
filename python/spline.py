import numpy as np

class CubicSpline:
    # x and y are numpy arrays
    def __init__(self, x, y):
        self.x = x
        self.y = y

        n = len(x)
        
        M = np.zeros((n, n))
        b = np.zeros(n)

        M[0, 0] = 2.0/(x[1] - x[0])
        M[0, 1] = 1.0/(x[1] - x[0])
        b[0] = 3.0*(y[1] - y[0])/(x[1] - x[0])**2

        for i in range(1, n-1):
            for j in range(i-1, i+2):
                if j == i-1:
                    M[i, j] = 1.0/(x[i] - x[i-1])
                elif j == i:
                    M[i, j] = 2.0*(1.0/(x[i] - x[i-1]) + 1.0/(x[i+1] - x[i]))
                elif j == i+1:
                    M[i, j] = 1.0/(x[i+1] - x[i])
            
            b[i] = 3.0*(y[i+1] - y[i])/(x[i+1] - x[i])**2 + 3.0*(y[i] - y[i-1])/(x[i] - x[i-1])**2
        
        M[n-1, n-2] = 1.0/(x[n-1] - x[n-2])
        M[n-1, n-1] = 2.0/(x[n-1] - x[n-2])
        b[n-1] = 3.0*(y[n-1] - y[n-2])/(x[n-1] - x[n-2])**2

        k = np.linalg.solve(M, b)

        c = np.zeros(n-1)
        d = np.zeros(n-1)

        for i in range(n-2):
            c[i] = k[i]*(x[i+1] - x[i]) - (y[i+1] - y[i])
            d[i] = -k[i+1]*(x[i+1] - x[i]) + (y[i+1] - y[i])

        self.k = k
        self.c = c
        self.d = d

    
    def __call__(self, x):
        i = np.searchsorted(self.x, x, "right") - 1

        t = (x - self.x[i])/(self.x[i+1] - self.x[i])
        return (1.0 - t)*self.y[i] + t*self.y[i+1] + t*(1.0 - t)*(self.c[i]*(1.0 - t) + self.d[i]*t)