import matplotlib
import matplotlib.pyplot as plt
import os
import sys
import numpy as np
import string
import math
import operator
from scipy import ndimage

from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm # color map




matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = 12

PATH_EXAMPLE = os.path.abspath(os.path.join(__file__, "../"))
PATH_THIS    = os.path.abspath(__file__)
PATH_PROJECT  = os.path.abspath(os.path.join(PATH_EXAMPLE, "../.."))
PATH_PYSCRIPT = os.path.abspath(os.path.join(PATH_PROJECT, "pyscript"))

sys.path.append(PATH_PYSCRIPT)
print PATH_PYSCRIPT
import read

def draw_contour():
    fig = plt.figure(figsize=(5, 4))
    ax  = fig.add_subplot(111, projection = '3d')

    """
    Set labels
    """
    plt.xlabel(r'x')
    plt.ylabel(r'y')


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
    cs   = read.CenterScalar("center_phi")
    assert(cs.get_dim() == 2)
    x, y = np.meshgrid(cs.get_coo_x(), cs.get_coo_y())
    z    = cs.get_mat_val()

    #fig = plt.contourf(x, y, z)
    p = ax.plot_surface(x, y, z, cmap=cm.coolwarm, antialiased = False)

    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    angle = -120
    ax.view_init(30, angle)
    #plt.clabel(fig, inline=1, fontsize=10, fmt="%.3f")
    # make a colorbar for the contour lines
    # CB = p.colorbar(fig, shrink=0.8, extend='both')

    #plt.grid(True)
    #plt.axes().set_aspect('equal', 'datalim')
    #plt.axes().set_aspect('equal')
    plt.tight_layout()
    plt.savefig("fig_phi.pdf")
    # plt.show()
    # rotate the axes and update


def draw_residual():
    plt.figure(figsize=(5, 4))

    """
    Set labels
    """
    plt.xlabel(r'iteration')
    plt.ylabel(r'Residual $||\mathbf{r}||/||\mathbf{b}||$')

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
    cs = read.read_residual_file("./residual")
    assert(len(cs) == 1)
    x = read.col(cs[0], 0)
    y = read.col(cs[0], 1)

    plt.plot(x, y, marker=".")
    plt.axes().set_yscale("log")
    #plt.axes().set_aspect('equal')
    plt.tight_layout()
    plt.savefig("fig_residual.pdf")

def main():
    draw_contour() 
    draw_residual()

if __name__ == '__main__':
    main()
