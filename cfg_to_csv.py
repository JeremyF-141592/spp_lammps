"""
Converts .cfg files to .csv, while merging the (mux, muy) properties describing
the orientation vector into a single variable "theta".

"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from glob import glob


def cfg_to_csv(path_to_file, destination_folder):
    export_path = os.path.join(destination_folder, path_to_file.split("/")[-1][:-4] + ".csv")
    if os.path.exists(export_path): return "File exists : " + export_path
    with open(path_to_file, "r") as f:
        big_list = list()
        register = False
        part = 0
        columns = list()
        frame = 0
        for line in f.readlines():
            if "ATOMS id x y" in line:
                register = True
                columns = line.split()[2:]
                part = 0
                continue
            elif "ITEM" in line and register is True:
                register = False
                frame += 1
                continue
            if register:
                row = [float(p) for p in line.split()]
                if "mux" in columns:
                    theta = np.arctan2(row[columns.index("muy")], row[columns.index("mux")])
                    remove_elements = [row[columns.index("muy")], row[columns.index("mux")]]
                    row.remove(remove_elements[0])
                    row.remove(remove_elements[1])
                    row.append(theta)
                part = row[columns.index("id")]
                big_list.append(row + [frame, part])
    if "mux" in columns:
        columns.remove("mux")
        columns.remove("muy")
        df = pd.DataFrame(big_list, columns=columns + ["theta", "frame", "particle"])
    else:
        df = pd.DataFrame(big_list, columns=columns + ["frame", "particle"])
    df.to_csv(export_path, index=False)
    return export_path


if __name__ == "__main__":

    for p in glob("./*.cfg"):
        p2 = cfg_to_csv(p, "./CSV/")
        print(p2)


