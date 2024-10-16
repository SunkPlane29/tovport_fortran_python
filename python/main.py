from constants import *
from spline import CubicSpline
from diff import solve_diffeq
from tov import solve_tov, solve_mrdiagram, initialize_eos

import time

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def main():
    eosfilename = "../eos_2.csv"

    initialize_eos(eosfilename)
    p0 = np.linspace(0.5*MEVFM3_TO_PRESSURE_UNIT, 600.0*MEVFM3_TO_PRESSURE_UNIT, 200)
    mrdiagram = solve_mrdiagram(p0, 12)

    #create dataframe without index and header
    df = pd.DataFrame(mrdiagram)
    df.to_csv("out/mrdiagram.dat", sep=' ', index=False, header=False)

if __name__ == '__main__':
    main()
