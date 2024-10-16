import matplotlib.pyplot as plt
import scienceplots
import pandas as pd

def main():
    eos_df = pd.read_csv("eos_1.csv", header=None, names=["p", "e"], skipinitialspace=True)

    #NOTE: I dont know how to read fortran data files
    fortran_df = pd.read_csv("fortran/out/mrdiagram.dat", sep=" ", header=None, names=["p0", "m", "r", "nan"], skipinitialspace=True)
    python_df = pd.read_csv("python/out/mrdiagram.dat", sep=" ", header=None, names=["p0", "m", "r"], skipinitialspace=True)
    julia_df = pd.read_csv("julia/out/mrdiagram.dat", sep=" ", header=None, names=["p0", "m", "r"], skipinitialspace=True)

    plt.style.use("science")

    fig, ax = plt.subplots()
    ax.plot(eos_df["p"], eos_df["e"], label="EOS", linestyle="-", color="red", linewidth=1)
    ax.set_xlabel("$P$ [MeV/fm$^3$]")
    ax.set_ylabel("$\epsilon$ [MeV/fm$^3$]")
    ax.legend()
    plt.savefig("out/eos.png", dpi=400)

    fig, ax = plt.subplots()
    # ax.plot(fortran_df["r"], fortran_df["m"], label="Fortran", linestyle="-", color="grey", linewidth=1)
    ax.plot(julia_df["r"], julia_df["m"], label="Julia", linestyle=(0, (5, 2)), color="green", linewidth=1)
    ax.plot(python_df["r"], python_df["m"], label="Python", linestyle=(0, (1, 3)), color="blue", linewidth=1)
    ax.set_xlabel("$R$ [km]")
    ax.set_ylabel("$M$ [M$_\odot$]")
    ax.legend()
    plt.savefig("out/mrdiagramjoint.png", dpi=400)

if __name__ == '__main__':
    main()
