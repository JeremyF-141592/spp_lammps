import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob


def w_distribution_animation(df, name="", step=10, pause=0.01):
    plt.show()
    for t in range(1, df["frame"].max(), step):
        frame = df[df["frame"] == t]
        plt.hist(0.5 * (1.+np.cos(frame["phi0"])), bins=20, density=True)
        #plt.hist(frame["phi0"], bins=20, density=True)
        plt.title(name + "\n" + str(t))
        plt.ylim(0, 4)
        plt.xlim(0, 1)
        plt.pause(pause)
        plt.clf()


def q_distribution_animation(df, name="", step=10, pause=0.01):
    plt.show()
    for t in range(1, df["frame"].max(), step):
        frame = df[df["frame"] == t]
        plt.hist(frame["qreward"], bins=np.linspace(0, 1, num=25))
        #plt.hist(frame["phi0"], bins=20, density=True)
        plt.title(name + "\n" + str(t))
        plt.ylim(1, 1000)
        plt.xlim(0, 1)
        plt.yscale('log')
        plt.pause(pause)
        plt.clf()


def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs."""
    return np.isnan(y), lambda z: z.nonzero()[0]


def get_mag(df, qtty="phi0", step=10):
    m0 = np.zeros(df["frame"].max() - 1) * np.nan
    for t in range(1, df["frame"].max(), step):
        frame = df[df["frame"] == t]
        m0[t - 1] = (np.sum(np.cos(frame[qtty])) ** 2 + np.sum(np.sin(frame[qtty])) ** 2) / (1000 ** 2)
    if step > 1:
        nans, x = nan_helper(m0)
        m0[nans] = np.interp(x(nans), x(~nans), m0[~nans])
    return m0


def learning_total(df, step=10):
    density = np.zeros(df["frame"].max() - 1) * np.nan
    up = df.loc[df["x"] < 25]
    tot = up.loc[up["x"] > -25]
    for t in range(1, df["frame"].max(), step):
        frame = tot[tot["frame"] == t]
        density[t - 1] = len(frame) / 1000
    if step > 1:
        nans, x = nan_helper(density)
        density[nans] = np.interp(x(nans), x(~nans), density[~nans])
    return density


def plot_m0(files, step=10):
    etas = list()
    for p in files:
        etas.append(float(p.split("/")[-1].split("_")[2]))
    files = [x for _, x in sorted(zip(etas, files))]
    for p in files:
        df = pd.read_csv(p)
        m0 = get_mag(df, "phi0", step=step)
        eta = p.split("/")[-1].split("_")[2]
        plt.plot(m0, label=eta)
    plt.ylim(0, 1)
    plt.legend()
    plt.title(r"$\frac{1}{N^2} (\Sigma (\cos \phi^0)^2 + \Sigma (\sin \phi^0)^2)$")
    plt.grid()
    plt.show()


def plot_learning_total(files, step=10):
    etas = list()
    alphaq = list()
    for p in files:
        etas.append(float(p.split("/")[-1].split("_")[2]))
        alphaq.append(float(p.split("/")[-1].split("_")[7]))
    files = [x for _, x in sorted(zip(alphaq, files))]
    indices = np.argsort(np.argsort(alphaq))
    for i in range(len(files)):
        df = pd.read_csv(files[i])
        tot = learning_total(df, step=step)
        color = (indices[i] / np.max(indices), 0, 1-indices[i] / np.max(indices))
        plt.plot(tot, label=f"{alphaq[i]} - {etas[i]}", color=color)
    plt.ylim(0, 1)
    plt.legend()
    plt.title("Total particles in center")
    plt.grid()
    plt.show()


if __name__ == "__main__":
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "Helvetica"
    })

    plot_learning_total(glob("/home/morphoswarm/Documents/spp_lammps/learning/simple_learning/evo2D_1000_0_*.csv"), step=100)
    df = pd.read_csv("/home/morphoswarm/Documents/spp_lammps/learning/simple_learning/evo2D_1000_0.01_0.001_1_0_1000_0.1_1_2_2349.csv")
    q_distribution_animation(df, "0.001 - 10e-2", step=50)
    plot_m0(glob("/home/morphoswarm/Documents/spp_lammps/learning/uniform_eta/*.csv"))
