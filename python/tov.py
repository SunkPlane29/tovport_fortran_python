from spline import CubicSpline
from constants import *
from diff import solve_diffeq

from multiprocessing import Pool

import pandas as pd
import numpy as np

def initialize_eos(eosfilename):
    global eos
    eos = EoS(eosfilename)

class EoS:
    def __init__(self, eosfilename):
        df = pd.read_csv(eosfilename, sep=',', header=None)
        P = np.array(df[0].values) * MEVFM3_TO_PRESSURE_UNIT
        e = np.array(df[1].values) * MEVFM3_TO_PRESSURE_UNIT

        self.cs = CubicSpline(P, e)

    
    def __call__(self, P):
        return self.cs(P)


def pressure_diffeq(r, x):
    P = x[0]
    M = x[1]

    return -(eos(P)*M/r**2)*(1 + P/eos(P))*(1 + 4*pi*r**3*P/M)*(1 - 2*M/r)**(-1)

def mass_diffeq(r, x):
    P = x[0]
    M = x[1]

    return 4*pi*r**2*eos(P)

def tov_diffeq(r, x): 
    return np.array([pressure_diffeq(r, x), mass_diffeq(r, x)])

def solve_tov(p0):
    t0 = 1.0e-8
    x0 = np.array([p0, 1.0e-24])
    h = 1.0*SI_TO_LENGTH_UNIT
    maxiter = 100000

    f = lambda t, x: tov_diffeq(t, x)
    terminate = lambda i, t, x: x[0] <= 0

    sol = solve_diffeq(f, t0, x0, h, maxiter, terminate)
    sol[:, 0] = sol[:, 0]*LENGTH_UNIT_TO_SI*1.0e-3
    sol[:, 1] = sol[:, 1]*PRESSURE_UNIT_TO_MEVFM3

    return sol

#returns [p0, mass, radius]
def custom_solve_tov(p0):
    sol = solve_tov(p0)
    return [p0, sol[-1, 2], sol[-1, 0]]

def solve_mrdiagram(p0, nthreads):
    pool = Pool(processes=nthreads)

    mrdiagram = np.array(pool.map(custom_solve_tov, p0)).reshape(len(p0), 3)
    pool.close()
    pool.join()
    return mrdiagram
    