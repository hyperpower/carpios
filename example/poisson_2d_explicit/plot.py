import matplotlib
import matplotlib.pyplot as plt
import os
import sys
import numpy as np
import string
import math
import operator
from scipy import interpolate
from scipy import ndimage

PATH_EXAMPLE  = os.path.abspath(os.path.join(__file__, "../"))
PATH_THIS     = os.path.abspath(__file__)
PATH_RESULT   = os.path.abspath("./result")
PATH_FIG      = os.path.abspath("./fig")
PATH_PROJECT  = os.path.abspath(os.path.join(PATH_EXAMPLE, "../.."))
PATH_PYSCRIPT = os.path.abspath(os.path.join(PATH_PROJECT, "pyscript"))

sys.path.append(PATH_PYSCRIPT)
import read


def _col(matrix, i):
    return [row[i] for row in matrix]

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = 12

def file_name(namedir, namevar):
    res = []
    files = [f for f in os.listdir(namedir) if os.path.isfile(os.path.join(namedir, f))]
    for f in files:
        spf = f.split("_")
        if spf[0] == namevar:
            res.append(spf)
    return res 

def plot_phi_one(strstep, strtime):
    fnp = PATH_RESULT + "/phi_" + strstep + "_" + strtime  

    plt.figure(figsize=(5, 5))
    
    x_st = 0.0
    x_ed = 1.0
    y_st = 0.0
    y_ed = 1.0

    plt.xlim([x_st, x_ed])
    plt.ylim([y_st, y_ed])

    plt.xlabel(r'x')
    plt.ylabel(r'y')

    pd   = read.PointData(fnp)
    arrx = pd.get_coo_x()
    arry = pd.get_coo_y()

    P    = pd.get_mat_val()

    X, Y = np.meshgrid(arrx, arry)

    le = np.linspace(0.0, 10.0, 10)

    Q = plt.contour(X, Y, P, levels = le, cmap="bwr")
    plt.colorbar()

    plt.text(0.0, 1.01, "Time = "+ "%.5f" % float(strtime))
    
    #plt.grid(True)
    plt.axes().set_aspect('equal')
    plt.tight_layout()
    name = PATH_FIG + "/phi_" + "%05d" % int(strstep) +".png"
    print name
    plt.savefig(name)
    plt.close()


def plot():
    matfu = file_name(PATH_RESULT, "phi")
    for one in matfu:
    #one = 
        print "Draw step : ", 
        print one[1] + " "
        plot_phi_one(one[1], one[2])

    #os.system("convert -delay 20 -loop 0 ./fig/phi_*.png phi.gif")

def main():
    plot() 

if __name__ == '__main__':
    main()
