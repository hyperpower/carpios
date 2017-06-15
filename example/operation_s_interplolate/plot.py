import matplotlib
import matplotlib.pyplot as plt
import os
import sys
import numpy as np
import string
import math
import operator
from scipy import ndimage

PATH_EXAMPLE  = os.path.abspath(os.path.join(__file__, "../"))
PATH_THIS     = os.path.abspath(__file__)
PATH_PROJECT  = os.path.abspath(os.path.join(PATH_EXAMPLE, "../.."))
PATH_PYSCRIPT = os.path.abspath(os.path.join(PATH_PROJECT, "pyscript"))

sys.path.append(PATH_PYSCRIPT)
print PATH_PYSCRIPT
import read


def _col(matrix, i):
    return [row[i] for row in matrix]

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = 12

PATH_EXAMPLE = os.path.abspath(os.path.join(__file__, "../"))
PATH_THIS    = os.path.abspath(__file__)

def draw():
    plt.figure(figsize=(4, 4))

    """
    Set labels
    """
    plt.xlabel(r'x')
    plt.ylabel(r'phi')

    """
    Set range
    """
    x_st = 3e-2
    x_ed = 100
    y_st = 1e-2
    y_ed = 20

    #plt.xlim([x_st, x_ed])
    #plt.ylim([y_st, y_ed])
    #plt.xscale('log')
    #plt.yscale('log')

    """
    Data part
    """
    # read file as point data
    pd  = read.PointData("center_phi")
    mat = pd._mat;
    plt.plot(_col(mat,0), _col(mat, 3), marker=".")

    res  = read.PointData("result")
    mat = res._mat
    lres, = plt.plot(_col(mat,0), _col(mat, 3), marker=".", linestyle="none")

    res  = read.PointData("exact")
    mat = res._mat
    lexa, = plt.plot(_col(mat,0), _col(mat, 3), marker=".", linestyle="none")

    plt.legend(
        [lres, lexa],
        ["Result", "Exact"]
        )

    plt.grid(True)
    #plt.axes().set_aspect('equal')
    plt.tight_layout()
    plt.savefig("fig_phi.pdf")
    # plt.show()


def main():
    draw() 

if __name__ == '__main__':
    main()
