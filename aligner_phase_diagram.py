import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
from matplotlib.colors import LogNorm


def order_parameter(df, timesteps=100):
    """Magnitude of mean velocity, last N timesteps"""
    mean_vx = np.mean(df["vx"].loc[df["frame"] > np.max(df["frame"]) - timesteps])
    mean_vy = np.mean(df["vy"].loc[df["frame"] > np.max(df["frame"]) - timesteps])
    return np.sqrt(mean_vy**2 + mean_vx**2)


def phase_diagram(files, density=0.01, timesteps=20):
    dic = dict()

    scattering_rate = 4*0.01 / np.pi
    # dic key is (noise, density), value is (order parameter,number of independent runs)
    xlabels = ["",]
    for p in files:
        param = order_parameter(pd.read_csv(p), timesteps=timesteps)
        print(param)
        noise = float(p.split("/")[-1].split("_")[2]) / scattering_rate
        alpha = float(p.split("/")[-1].split("_")[-2]) * 10
        seed = int(p.split("/")[-1].split("_")[-1][:-4])

        if (noise, alpha) in dic.keys():
            independant_runs = dic[(noise, alpha)][1] + 1
            new_value = dic[(noise, alpha)][0] * dic[(noise, alpha)][1] + param
            new_value /= independant_runs
            dic[(noise, alpha)] = (new_value, independant_runs)
        else:
            dic[(noise, alpha)] = (param, 1)

    x = list()
    y = list()
    z = list()

    for key in dic.keys():
        y.append(key[0])
        x.append((key[1]))
        z.append(dic[key][0])
    print(x, y, z)

    plt.scatter(x, y, s=150, c=z, cmap='seismic', marker="d")
    for i in range(len(x)):
        plt.annotate(str(x[i]), (x[i], y[i]+0.05))
    ax = plt.gca()
    plt.xlabel("alpha", fontsize=18)
    plt.ylabel("D \ lambda", fontsize=18)
    plt.axhline(0, color="k")

    plt.ylim(-0.05, 0.31)
    plt.xlim(0.01, 100)
    ax.set_xscale('log')
    # ax.set_facecolor((0.8, 0.8, 0.8))
    plt.colorbar()
    plt.draw()

    # ax.set_xticklabels(xlabels)

    plt.show()


def annealing_plot(df, D):
    v = df[["frame", "vx", "vy"]]
    v = v.sort_values("frame")
    parameters = np.zeros(3600)
    for i in range(3600):
        vframe = v[v["frame"] == i]
        parameters[i] = np.sqrt(np.mean(vframe["vx"])**2 + np.mean(vframe["vy"])**2)
        print(i)
    plt.plot(parameters)
    plt.show()
    x_forward = np.linspace(0, D, num=1200)
    x_backward = np.linspace(D, 0, num=1200)
    plt.plot(x_forward, parameters[1200:2400], "r")
    plt.plot(x_backward, parameters[2400:], "b")
    plt.xlabel("D / lambda")
    plt.ylabel("phi")
    plt.show()


if __name__ == "__main__":
    annealing_plot(pd.read_csv("/home/morphoswarm/Documents/spp_lammps/examples/annealing/alpha1/spp_1000_annealing0.00635_0.1_0.1_1236.csv"), 0.5)
    annealing_plot(pd.read_csv("/home/morphoswarm/Documents/spp_lammps/examples/annealing/alpha01/spp_1000_annealing0.00127_0.1_0.01_1236.csv"), 0.1)
    phase_diagram(glob("/home/morphoswarm/Documents/spp_lammps/examples/aligner_collective/*.csv"))
