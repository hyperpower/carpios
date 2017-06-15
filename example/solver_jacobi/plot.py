import matplotlib
import matplotlib.pyplot as plt
import os
import numpy as np
import string
import math
import operator
from scipy import ndimage


def read_residual(fn):
    file = open(fn, "r") 
    content = file.readlines()
    arr = []
    for row in content[2:]:
        row.strip()
        num = float(row)
        arr.append(num)
    return arr

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
    plt.ylabel(r'Residual')
    plt.xlabel(r'Number of iteration')

    """
    Set range
    """
    x_st = 3e-2
    x_ed = 100
    y_st = 1e-2
    y_ed = 20

    plt.xlim([1, 1e5])
    #plt.ylim([y_st, y_ed])
    plt.xscale('log')
    plt.yscale('log')

    """
    Data part
    """
    omega = [1.0, 1.2, 1.5, 1.7, 1.9]
    lname = ["Jacobi"]
    for o in omega:
        n = "SOR_" + "%.2f" % o
        lname.append(n)
    mat = []
    for n in lname:
        arr = read_residual(n + "_residual.txt")
        mat.append(arr)
    
    #assert(int(d) == 1)
    lleg = []
    for row in mat:
        lidx = []
        for i in range(0, len(row)):
            lidx.append(i) 
        
        leg, = plt.plot(lidx, row)
        lleg.append(leg)


    plt.legend(lleg, [ (n.replace("_", " $\omega$ = ") ) for n in lname])

    plt.grid(True)
    #plt.axes().set_aspect('equal')
    plt.tight_layout()
    plt.savefig("residual_comp.pdf")
    # plt.show()


def main():
    draw() 

if __name__ == '__main__':
    main()
