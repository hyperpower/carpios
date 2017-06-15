import matplotlib
import matplotlib.pyplot as plt
import os
import numpy as np
import string
import math
import operator
from scipy import ndimage


def read_center_scalar(fn):
    file = open(fn, "r") 
    content = file.readlines()
    mat = []
    for row in content:
        nr = [x.strip() for x in row.split(',')]
        mat.append(nr)
    fl   = mat[0]
    sfl  = fl[0].split(' ')
    nr   = (sfl[1].split(':'))[1]
    dim  = (sfl[2].split(':'))[1]
    nx   = (sfl[3].split(':'))[1]
    ny   = (sfl[4].split(':'))[1]
    nz   = (sfl[5].split(':'))[1]
    mat  = mat[1:]
    return mat, nr, dim, nx, ny, nz

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
    mat, nr, d ,nx, ny, nz = read_center_scalar("center_phi")
    assert(int(d) == 1)
    plt.plot(_col(mat,0), _col(mat, 3), marker=".")

    plt.grid(True)
    #plt.axes().set_aspect('equal')
    plt.tight_layout()
    plt.savefig("fig_phi.pdf")
    # plt.show()


def main():
    draw() 

if __name__ == '__main__':
    main()
