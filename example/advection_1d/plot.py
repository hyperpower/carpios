import matplotlib
import matplotlib.pyplot as plt
import os, sys
import numpy as np
import string
import math
import operator
from scipy import ndimage
import multiprocessing

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

def split(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out

def plot_one(strstep, strtime):
    print "Draw : ", strstep, " ", strtime

    plt.figure(figsize=(6, 4))

    """
    Set labels
    """
    plt.xlabel(r'x')
    plt.ylabel(r'phi')

    """
    Set range
    """
    x_st = 0
    x_ed = 200
    y_st = -0.5
    y_ed = 1.5

    plt.xlim([x_st, x_ed])
    plt.ylim([y_st, y_ed])
    #plt.xscale('log')
    #plt.yscale('log')

    """
    Data part
    """
    fne = PATH_RESULT + "/exact_" + strstep + "_" + strtime  

    scheme = ["upwind1", "center", "center4"]

    fnv = []
    for s in scheme:
        fnv.append(PATH_RESULT + "/phi_" + s + "_" + strstep + "_" + strtime)  

    pe   = read.PointData(fne)
    arrx = pe.get_coo_x()
    arre  = pe.get_arr_val()
    arrv = []
    for f in fnv: 
        pv   = read.PointData(f)
        arrv.append(pv.get_arr_val())
    
    plt.plot(arrx, arre) 
    llg = []
    for arr in arrv:
        lg, = plt.plot(arrx, arr, marker=".") 
        llg.append(lg)

    plt.text(10, 1.25, "Time = "+ "%.2f" % float(strtime))
    plt.text(10, 1.00, "Step = "+ "%04d" % float(strstep))

    plt.legend(llg, scheme, loc= 'upper right')

    plt.grid(True)
    #plt.axes().set_aspect('equal')
    plt.tight_layout()
    plt.savefig(PATH_FIG + "/comp_" + "%06d" % int(strstep) +".png")
    plt.close()
    # plt.show()


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

def plot_error1():
    plt.figure(figsize=(4, 4))

    """
    Set labels
    """
    plt.ylabel(r'Error$_1$')
    plt.xlabel(r'Step')

    """
    Set range
    """
    x_st = 0
    x_ed = 200
    y_st = -0.5
    y_ed = 1.5

    #plt.xlim([x_st, x_ed])
    #plt.ylim([y_st, y_ed])
    #plt.xscale('log')
    #plt.yscale('log')

    """
    Data part
    """
    scheme = ["upwind1", "center", "center4"]
    mat = []
    for s in scheme:
        m = read_error_file(s)
        mat.append(m)
    llg = []
    for i, v in enumerate(scheme):
        arrx  = _col(mat[i], 0)
        arre1 = _col(mat[i], 2)
        arre2 = _col(mat[i], 3)
        arrei = _col(mat[i], 4)

        lg, = plt.plot(arrx, arre1) 
        llg.append(lg)

    #plt.text(10, 1.25, "Time = "+ "%.2f" % float(strtime))
    #plt.text(10, 1.00, "Step = "+ "%04d" % float(strstep))

    plt.legend(llg, scheme, loc= 'upper left')

    plt.grid(True)
    #plt.axes().set_aspect('equal')
    plt.tight_layout()
    plt.savefig("error_1.png")
    plt.close()
    # plt.show()

def plot_error2():
    plt.figure(figsize=(4, 4))

    """
    Set labels
    """
    plt.ylabel(r'Error$_2$')
    plt.xlabel(r'Step')

    """
    Set range
    """
    x_st = 0
    x_ed = 200
    y_st = -0.5
    y_ed = 1.5

    #plt.xlim([x_st, x_ed])
    #plt.ylim([y_st, y_ed])
    #plt.xscale('log')
    #plt.yscale('log')

    """
    Data part
    """
    scheme = ["upwind1", "center", "center4"]
    mat = []
    for s in scheme:
        m = read_error_file(s)
        mat.append(m)
    llg = []
    for i, v in enumerate(scheme):
        arrx  = _col(mat[i], 0)
        arre1 = _col(mat[i], 2)
        arre2 = _col(mat[i], 3)
        arrei = _col(mat[i], 4)

        lg, = plt.plot(arrx, arre2) 
        llg.append(lg)

    #plt.text(10, 1.25, "Time = "+ "%.2f" % float(strtime))
    #plt.text(10, 1.00, "Step = "+ "%04d" % float(strstep))

    plt.legend(llg, scheme, loc= 'upper left')

    plt.grid(True)
    #plt.axes().set_aspect('equal')
    plt.tight_layout()
    plt.savefig("error_2.png")
    plt.close()

def plot_errori():
    plt.figure(figsize=(4, 4))

    """
    Set labels
    """
    plt.ylabel(r'Error$_{\inf}$')
    plt.xlabel(r'Step')

    """
    Set range
    """
    x_st = 0
    x_ed = 200
    y_st = -0.5
    y_ed = 1.5

    #plt.xlim([x_st, x_ed])
    #plt.ylim([y_st, y_ed])
    #plt.xscale('log')
    #plt.yscale('log')

    """
    Data part
    """
    scheme = ["upwind1", "center", "center4"]
    mat = []
    for s in scheme:
        m = read_error_file(s)
        mat.append(m)
    llg = []
    for i, v in enumerate(scheme):
        arrx  = _col(mat[i], 0)
        arre1 = _col(mat[i], 2)
        arre2 = _col(mat[i], 3)
        arrei = _col(mat[i], 4)

        lg, = plt.plot(arrx, arrei) 
        llg.append(lg)

    #plt.text(10, 1.25, "Time = "+ "%.2f" % float(strtime))
    #plt.text(10, 1.00, "Step = "+ "%04d" % float(strstep))

    plt.legend(llg, scheme, loc= 'upper left')

    plt.grid(True)
    #plt.axes().set_aspect('equal')
    plt.tight_layout()
    plt.savefig("error_i.png")
    plt.close()

def plot_all():
    matfu = file_name(PATH_RESULT, "exact")
    print len(matfu)
    matfc = []
    for one in matfu:
        matfc.append(one)

    #multiprocessing.freeze_support()
    pool = multiprocessing.Pool()
    cpus = multiprocessing.cpu_count() / 2
    results = []
    cmatfs = split(matfc, cpus)

    print "Thread num : ", len(cmatfs)
    for i in xrange(0, cpus):
        mat = cmatfs[i]
        for one in mat:
            result = pool.apply_async(plot_one, args=(one[1], one[2],))
            results.append(result)
    
    pool.close()
    pool.join()
    
    os.system("convert -delay 5 -loop 0 ./fig/comp_*.png comp.gif")


def main():
    #stri = "480"
    #strt = "48"
    #plot_one(stri, strt) 
    
    plot_all()
    plot_error1()
    plot_error2()
    plot_errori()

if __name__ == '__main__':
    main()


