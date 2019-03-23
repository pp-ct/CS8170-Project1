import numpy as np
from os.path import dirname, abspath
import sys
DIRECTORY = dirname(dirname(abspath(__file__)))
sys.path.append(DIRECTORY)

from lib import library


DATA_DIR = "../data/"

def distance_matrix(filename, target_name):
    input = np.genfromtxt(filename)

    coord = input[:, :2]
    data = input[:, 2:]

    x = input[:, :1].flatten()
    y = input[:, 1:2].flatten()
    std = np.mean(input[:, 3:5], axis=1).flatten()
    mean = input[:, 2].flatten()

    lengthX = int(np.max(x))
    lengthY = int(np.max(y))
    length = max(lengthX, lengthY)


    mat = np.full((length, length, 2), np.nan)
    for i in range(len(x)):
        x_coord = int(x[i] - 1)
        y_coord = int(y[i] - 1)

        mat[x_coord][y_coord][0] = mean[i]
        mat[x_coord][y_coord][1] = std[i]

        mat[y_coord][x_coord][0] = mean[i]
        mat[y_coord][x_coord][1] = std[i]

    library.output_distance_matrix(DATA_DIR, mat[:, :, 0], prefix=target_name)
    return mat



