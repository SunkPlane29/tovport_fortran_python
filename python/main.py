from constants import *
from spline import CubicSpline
from diff import solve_diffeq
from tov import solve_tov, solve_mrdiagram, initialize_eos

import time

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def main():
    eosfilename = "../eos_1.csv"

    initialize_eos(eosfilename)
    p0i = 5.0
    p0f = 800.0
    n = 200
    p0log = np.linspace(np.log(p0i), np.log(p0f), n)
    p0 = np.exp(p0log)*MEVFM3_TO_PRESSURE_UNIT
    mrdiagram = solve_mrdiagram(p0, 12)

    #create dataframe without index and header
    df = pd.DataFrame(mrdiagram)
    df.to_csv("out/mrdiagram.dat", sep=' ', index=False, header=False)

if __name__ == '__main__':
    main()
