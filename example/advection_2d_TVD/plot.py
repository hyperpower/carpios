import matplotlib
import matplotlib.pyplot as plt
import os, sys
import numpy as np
import string
import math
import operator
from scipy import ndimage
from scipy import interpolate

import threading
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm # color map

from run import scheme as SCHEME
import copy



def _col(matrix, i):
    return [row[i] for row in matrix]

def chunks(l, n):
    n = max(1, n)
    return [l[i:i+n] for i in xrange(0, len(l), n)]

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = 12

PATH_EXAMPLE  = os.path.abspath(os.path.join(__file__, "../"))
PATH_THIS     = os.path.abspath(__file__)
PATH_RESULT   = os.path.abspath("./result")
PATH_FIG      = os.path.abspath("./fig")
PATH_PROJECT  = os.path.abspath(os.path.join(PATH_EXAMPLE, "../.."))
PATH_PYSCRIPT = os.path.abspath(os.path.join(PATH_PROJECT, "pyscript"))

sys.path.append(PATH_PYSCRIPT)
import read

def file_name(namedir, namevar):
    res = []
    files = [f for f in os.listdir(namedir) if os.path.isfile(os.path.join(namedir, f))]
    for f in files:
        spf = f.split("_")
        if spf[0] == namevar:
            res.append(spf)
    return res 

def plot_one(scheme, strstep, strtime):

    fig = plt.figure(figsize=(6, 4))
    #ax  = fig.add_subplot(111, projection = '3d')

    """
    Set labels
    """
    plt.xlabel(r'x')
    plt.ylabel(r'y')

    """
    Set range
    """
    x_st = 0
    x_ed = 1
    y_st = 0
    y_ed = 1

    plt.xlim([x_st, x_ed])
    plt.ylim([y_st, y_ed])
    #plt.zlim([-0.3, 1.3])
    #plt.xscale('log')
    #plt.yscale('log')

    """
    Data part
    """
    fne = PATH_RESULT + "/exact_" + strstep + "_" + strtime  

    #scheme = ["upwind2", "VanLeer", "superbee", "WAHYD"]

    fnv = PATH_RESULT + "/phi_" + scheme + "_" + strstep + "_" + strtime

    pe   = read.PointData(fne)
    arrx = pe.get_coo_x()
    arry = pe.get_coo_y()
    X, Y = np.meshgrid(arrx, arry)
    mate = pe.get_mat_val()
    pv  = read.PointData(fnv)
    mat = pv.get_mat_val()
    
    lg = plt.contour(X, Y, mat, cmap=cm.coolwarm, antialiased = False)
    plt.colorbar()

    fp = interpolate.interp2d(arrx, arry, mat, kind='linear')
    y1new = np.arange(0, 1.0, 1e-2)
    interp1  = fp(0.5, y1new)
    interp2  = fp(y1new, 0.5)

    plt.plot(interp1 * 0.5 + 0.5, y1new, color="g") 
    plt.plot(y1new, interp2 * 0.5 + 0.5, color="g") 


    #plt.text(10, 1.25, "Time = "+ "%.2f" % float(strtime))
    #plt.text(10, 1.00, "Step = "+ "%04d" % float(strstep))

    #plt.legend(llg, scheme, loc= 'upper right')
    #plt.grid(True)
    plt.axes().set_aspect('equal')
    plt.tight_layout()
    plt.savefig(PATH_FIG + "/"+ scheme +"_" + "%06d" % int(strstep) +".png")
    plt.close()
    #plt.show()

def read_error_file(scheme):
    filename = PATH_RESULT + "/error_" + scheme
    f = open(filename, 'r')
    content=[]
    for i, line in enumerate(f):
        line = line.strip()
        line = line.split()
        step = float(line[1])
        t    = float(line[3])
        e1   = float(line[5])
        e2   = float(line[7])
        ei   = float(line[9])
        content.append([step, t, e1, e2, ei])
    return content

def plot_all(scheme):
    matfu = file_name(PATH_RESULT, "phi")
    for one in matfu:
        print "Draw step : ",
        print one[2]
        plot_one(scheme, one[2], one[3])

    os.system("convert -delay 5 -loop 0 ./fig/" + scheme +"_*.png " + scheme + ".gif")

def plot_hline(strstep, strtime):
    fig = plt.figure(figsize=(4, 4))
    #ax  = fig.add_subplot(111, projection = '3d')

    """
    Set labels
    """
    plt.xlabel(r'x')
    plt.ylabel(r'phi')

    """
    Set range
    """
    x_st = 0
    x_ed = 1
    y_st = -0.1
    y_ed = 1.1

    plt.xlim([x_st, x_ed])
    plt.ylim([y_st, y_ed])
    #plt.zlim([-0.3, 1.3])
    #plt.xscale('log')
    #plt.yscale('log')

    """
    Data part
    """
    fne = PATH_RESULT + "/exact_" + strstep + "_" + strtime  
    fnv = []
    scheme = copy.copy(SCHEME)
    #scheme = ["upwind2", "VanLeer", "superbee", "WAHYD"]
    for s in scheme:
        fnv.append(PATH_RESULT + "/phi_" + s + "_" + strstep + "_" + strtime)

    pe   = read.PointData(fne)
    arrx = pe.get_coo_x()
    arry = pe.get_coo_y()
    X, Y = np.meshgrid(arrx, arry)
    mate = pe.get_mat_val()
    llg  = []

    for f in fnv:
        pv   = read.PointData(f)
        mat  = pv.get_mat_val()
        fp = interpolate.interp2d(arrx, arry, mat, kind='linear')
        y1new = np.arange(0, 1.0, 1e-2)
        interp1  = fp(0.5, y1new)
        interp2  = fp(y1new, 0.5)
        lg, = plt.plot(y1new, interp2) 
        llg.append(lg)

    fp = interpolate.interp2d(arrx, arry, mate, kind='linear')
    y1new = np.arange(0, 1.0, 1e-2)
    interp1  = fp(0.5, y1new)
    interp2  = fp(y1new, 0.5)
    lg, = plt.plot(y1new, interp2, color= "k") 
    llg.append(lg)

    scheme.append("Exact")
    plt.legend(llg, scheme, loc= 'upper right')
    #plt.grid(True)
    #plt.axes().set_aspect('equal')
    plt.tight_layout()
    plt.savefig("comp_h.png")
    plt.close()
    #plt.show()

def plot_vline(strstep, strtime):
    fig = plt.figure(figsize=(4, 4))
    #ax  = fig.add_subplot(111, projection = '3d')

    """
    Set labels
    """
    plt.xlabel(r'x')
    plt.ylabel(r'phi')

    """
    Set range
    """
    x_st = 0
    x_ed = 1
    y_st = -0.1
    y_ed = 1.1

    plt.xlim([x_st, x_ed])
    plt.ylim([y_st, y_ed])
    #plt.zlim([-0.3, 1.3])
    #plt.xscale('log')
    #plt.yscale('log')

    """
    Data part
    """
    fne = PATH_RESULT + "/exact_" + strstep + "_" + strtime  
    fnv = []
    scheme = copy.copy(SCHEME)
    #scheme = ["upwind2", "VanLeer", "superbee", "WAHYD"]
    for s in scheme:
        fnv.append(PATH_RESULT + "/phi_" + s + "_" + strstep + "_" + strtime)

    pe   = read.PointData(fne)
    arrx = pe.get_coo_x()
    arry = pe.get_coo_y()
    X, Y = np.meshgrid(arrx, arry)
    mate = pe.get_mat_val()
    llg  = []

    for f in fnv:
        pv   = read.PointData(f)
        mat  = pv.get_mat_val()
        fp = interpolate.interp2d(arrx, arry, mat, kind='linear')
        y1new = np.arange(0, 1.0, 1e-2)
        interp1  = fp(0.5, y1new)
        interp2  = fp(y1new, 0.5)
        lg, = plt.plot(y1new, interp1) 
        llg.append(lg)

    fp = interpolate.interp2d(arrx, arry, mate, kind='linear')
    y1new = np.arange(0, 1.0, 1e-2)
    interp1  = fp(0.5, y1new)
    interp2  = fp(y1new, 0.5)
    lg, = lg, = plt.plot(y1new, interp1, color= "k") 
    llg.append(lg)

    scheme.append("Exact")
    plt.legend(llg, scheme, loc= 'upper left')
    #plt.grid(True)
    #plt.axes().set_aspect('equal')
    plt.tight_layout()
    plt.savefig("comp_v.png")
    plt.close()


def plot_last():
    scheme = ["upwind2", "VanLeer", "superbee", "WAHYD"]
    matfu = file_name(PATH_RESULT, "phi")
    lstep = []
    for one in matfu:
        if(one[1] == scheme[1]):
            lstep.append([int(one[2]), one[0], one[1], one[2], one[3]])

    lstep = sorted(lstep,key=lambda x: (x[0]))
    
    last = lstep[-1]
    plot_hline(last[3], last[4])
    plot_vline(last[3], last[4])

    
    scheme = copy.copy(SCHEME)
    for s in scheme:
        plot_one(s, last[3], last[4])

def main():
    #tri = "970"
    #strt = "1.94"
    #plot_one("upwind2", stri, strt) 
    #scheme = ["upwind2", "VanLeer", "superbee", "WAHYD"]
    #plot_error1(scheme)
    #plot_error2(scheme)
    #plot_errori(scheme)
    #plot_all("upwind2")
    plot_last()

if __name__ == '__main__':
    main()


