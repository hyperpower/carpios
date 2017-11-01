import matplotlib
import matplotlib.pyplot as plt
import os
import sys
import numpy as np
import string
import math
import operator

PATH_EXAMPLE  = os.path.abspath(os.path.join(__file__, "../../"))
PATH_THIS     = os.path.abspath(__file__)
PATH_RESULT   = os.path.abspath("./result")
PATH_FIG      = os.path.abspath("./fig")
PATH_PROJECT  = os.path.abspath(os.path.join(PATH_EXAMPLE, "../.."))
PATH_PYSCRIPT = os.path.abspath(os.path.join(PATH_PROJECT, "pyscript"))

sys.path.append(PATH_PYSCRIPT)
import read

def draw_contour():
    plt.figure(figsize=(5, 4))

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
    cs = read.CenterScalar("center_phi")
    assert(cs.get_dim() == 2)
    x, y = np.meshgrid(cs.get_coo_x(), cs.get_coo_y())
    z    = cs.get_mat_val()

    fig = plt.contourf(x, y, z)
    plt.clabel(fig, inline=1, fontsize=10, fmt="%.3f")
    # make a colorbar for the contour lines
    CB = plt.colorbar(fig, shrink=0.8, extend='both')

    #plt.grid(True)
    plt.axes().set_aspect('equal', 'datalim')
    #plt.axes().set_aspect('equal')
    plt.tight_layout()
    fn = "fig_phi.pdf"
    plt.savefig(fn)
    print("Output fig : ", fn)
    # plt.show()

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
    fn = "fig_residual.pdf"
    plt.savefig(fn)
    print("Output fig : ", fn)

def main():
    draw_contour() 
    draw_residual()

if __name__ == '__main__':
    main()
