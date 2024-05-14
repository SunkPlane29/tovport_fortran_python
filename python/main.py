from constants import *
from spline import CubicSpline
from diff import solve_diffeq
from tov import solve_tov, solve_mrdiagram, initialize_eos

import time

import numpy as np
import matplotlib.pyplot as plt

def main():
    eosfilename = "../eos.csv"

    initialize_eos(eosfilename)
    p0 = np.linspace(1.0*MEVFM3_TO_PRESSURE_UNIT, 600.0*MEVFM3_TO_PRESSURE_UNIT, 200)
    start = time.time()
    mrdiagram = solve_mrdiagram(p0, 12)
    end = time.time()

    print("Elapsed time: ", end - start, " seconds")
    
    plt.plot(mrdiagram[:, 2], mrdiagram[:, 1])
    plt.xlabel("Radius (km)")
    plt.ylabel("Mass (solar mass)")
    plt.show()

if __name__ == '__main__':
    main()