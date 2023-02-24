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
    up = df.loc[df["x"] > 0]
    for t in range(1, df["frame"].max(), step):
        frame = up[up["frame"] == t]
        density[t - 1] = len(frame) / 1000
    if step > 1:
        nans, x = nan_helper(density)
        density[nans] = np.interp(x(nans), x(~nans), density[~nans])
    return density


def u_mean_std(df, step=10, intensity=0.0):
    res = np.zeros(df["frame"].max() - 1) * np.nan
    res2 = np.zeros(df["frame"].max() - 1) * np.nan
    for t in range(1, df["frame"].max(), step):
        frame = df[df["frame"] == t]
        w0 = 0.5 * (1 + np.cos(frame["phi0"]))
        w1 = 0.5 * (1 + np.cos(frame["phi1"]))
        u = 0.5 * (1 + np.tanh((w0 - intensity) / w1))
        res[t - 1] = np.mean(u)
        res2[t - 1] = np.std(u)
    if step > 1:
        nans, x = nan_helper(res)
        res[nans] = np.interp(x(nans), x(~nans), res[~nans])
        nans, x = nan_helper(res2)
        res2[nans] = np.interp(x(nans), x(~nans), res2[~nans])
    return res, res2


def phi_01(df, step=10):
    res = np.zeros(df["frame"].max() - 1) * np.nan
    res2 = np.zeros(df["frame"].max() - 1) * np.nan
    for t in range(1, df["frame"].max(), step):
        frame = df[df["frame"] == t]
        res[t - 1] = np.mean(frame["phi0"])
        res2[t - 1] = np.mean(frame["phi1"])
    if step > 1:
        nans, x = nan_helper(res)
        res[nans] = np.interp(x(nans), x(~nans), res[~nans])
        nans, x = nan_helper(res2)
        res2[nans] = np.interp(x(nans), x(~nans), res2[~nans])
    return res, res2


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
        print(p)
    plt.ylim(0, 1)
    plt.legend()
    plt.title(r"$\frac{1}{N^2} (\Sigma (\cos \phi^0)^2 + \Sigma (\sin \phi^0)^2)$")
    plt.grid()
    plt.show()


def phi0_distrib(files, bins=100, I=0.75):
    etas = list()
    for p in files:
        etas.append(float(p.split("/")[-1].split("_")[2]))
    files = [x for _, x in sorted(zip(etas, files))]
    for p, et in zip(files, etas):
        df = pd.read_csv(p)
        last_frame = df.loc[df["frame"] == df["frame"].max()]
        # w0 = 0.5*(1+np.cos(last_frame["phi0"]))
        # w1 = 0.5*(1+np.cos(last_frame["phi1"]))
        # w2 = 0.5*(1+np.cos(last_frame["phi2"]))
        # th = np.tanh((w0-I)/w1)
        # u = 0.5*(1+th)*w2 + 0.5*(1-th)*(1-w2)
        vmag = np.sqrt(last_frame["vx"]**2 + last_frame["vy"]**2)
        hist, bin_edges = np.histogram(vmag, bins=np.linspace(0, 1, num=bins), density=True)
        bin_center = (0.5 * (bin_edges + np.roll(bin_edges, -1)))[:-1]
        plt.plot(bin_center, hist, color="r" if et == 0.001 else "b", alpha=0.3)
        print(p)
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
        plt.plot(tot, label=f"{alphaq[i]} - {etas[i]}")
        print(f" {i} / {len(files)}")
    plt.axhline(0.5, color="k", linestyle="--")
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
    phi0_distrib(glob("/media/morphoswarm/MSSimulations/spp_lammps/learning/uniform_eta/CSV/**/*.csv", recursive=True))
    plot_m0(glob("/media/morphoswarm/MSSimulations/spp_lammps/learning/uniform_eta/CSV/**/*.csv", recursive=True), step=50)

    # df = pd.read_csv("/media/morphoswarm/MSSimulations/spp_lammps/learning/uniform_eta/evo2D_1000_0.001_0.001_1_0_1000_1000_1_2_2348.csv")
    #
    # xabs  = np.arange(df["frame"].max()-1)
    #
    # u_mean, u_std = u_mean_std(df)
    # plt.plot(xabs, u_mean, label="u", color="r")
    # plt.fill_between(xabs, u_mean - u_std, u_mean + u_std, color="r", alpha=0.1)
    # phi0, phi1 = phi_01(df)
    # plt.plot(xabs, phi0, label="phi0")
    # plt.plot(xabs, phi1, label="phi1")
    # plt.legend()
    # plt.grid()
    # plt.show()
    # alpha_dict = dict()
    # alpha_eta_dict = dict()
    # for p in glob("/media/morphoswarm/MSSimulations/spp_lammps/learning/learning_outlight/evo*.csv"):
    #     alpha = p.split("/")[-1].split("_")[6]
    #     eta = p.split("/")[-1].split("_")[2]
    #     if alpha not in alpha_dict.keys():
    #         alpha_dict[alpha] = list()
    #     if (alpha, eta) not in alpha_eta_dict.keys():
    #         alpha_eta_dict[(alpha, eta)] = list()
    #
    #     alpha_dict[alpha].append(p)
    #     alpha_eta_dict[(alpha, eta)].append(p)
    # plot_learning_total(alpha_dict["0.001"])
    # plot_learning_total(alpha_eta_dict[("0.01", "0")])  # bimodal
    # plot_learning_total(alpha_eta_dict[("0.01", "0.001")])  # bimodal
    # plot_learning_total(alpha_eta_dict[("0.01", "0.01")])  # bimodal
    # plot_learning_total(alpha_dict["0.01"])  # bimodal
    # plot_learning_total(alpha_dict["0.1"])
    # plot_learning_total(alpha_dict["1"])
    # plot_learning_total(alpha_dict["10"])
    #
    # df = pd.read_csv("/home/morphoswarm/Documents/spp_lammps/learning/simple_learning/evo2D_1000_0.01_0.001_1_0_1000_0.1_1_2_2349.csv")
    # q_distribution_animation(df, "0.001 - 10e-2", step=50)
    # plot_m0(glob("/home/morphoswarm/Documents/spp_lammps/learning/uniform_eta/*.csv"))
