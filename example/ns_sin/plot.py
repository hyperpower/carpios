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
import multiprocessing

PATH_EXAMPLE  = os.path.abspath(os.path.join(__file__, "../.."))
PATH_THIS     = os.path.abspath(__file__)
PATH_RESULT   = os.path.abspath("./result")
PATH_FIG      = os.path.abspath("./fig")
PATH_PROJECT  = os.path.abspath(os.path.join(PATH_THIS, "../../.."))
PATH_PYSCRIPT = os.path.abspath(os.path.join(PATH_PROJECT, "pyscript"))
PATH_PY       = os.path.abspath(os.path.join(PATH_EXAMPLE, "py"))

sys.path.append(PATH_PYSCRIPT)
import read
print PATH_PROJECT
sys.path.append(PATH_PY)
import botella_o as Botella

import run


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

def split(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out

def read_error_file(n):
    filename = PATH_RESULT + "/err" + n
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

def exact_u(x, y, t):
    return -math.cos(x) * math.sin(y) * math.exp(-2.0 * t)

def exact_v(x, y, t):
    return math.cos(y) * math.sin(x) * math.exp(-2.0 * t)

def exact_p(x, y, t):
    return -0.25* (math.cos(2 * x) + math.cos(2 * y)) * math.exp( - 4.0 * t)

def plot_velocity_exact(t):

    plt.figure(figsize=(5, 5))
    

    x_st = -math.pi / 2.0 
    x_ed = math.pi / 2.0
    y_st = -math.pi / 2.0 
    y_ed = math.pi / 2.0 

    plt.xlim([x_st, x_ed])
    plt.ylim([y_st, y_ed])

    plt.xlabel(r'x')
    plt.ylabel(r'y')

    arrx = np.linspace(x_st, x_ed, 40)
    arry = np.linspace(y_st, y_ed, 40)

    X, Y = np.meshgrid(arrx, arry)
    
    U = []
    for y in arry:
        row = []
        for x in arrx:
            v = exact_u(x, y, t)
            row.append(v)
        U.append(row)

    V = []
    for y in arry:
        row = []
        for x in arrx:
            v = exact_v(x, y, t)
            row.append(v)
        V.append(row)

    M = np.hypot(U, V)

    Q = plt.quiver(X, Y, U, V, M, cmap="bwr")
    plt.colorbar()

    #fu = interpolate.interp2d(arrx, arry, U, kind='linear')
    #fv = interpolate.interp2d(arrx, arry, V, kind='linear')

    #y1new = np.arange(-0.5, 0.5, 1e-2)

    #interu  = fu(0.0, y1new)
    #interv  = fv(y1new, 0.0)

    #plt.plot(interu * 0.5, y1new, color="g") 
    #plt.plot(y1new, interv * 0.5, color="g") 


    #plt.text(-0.5, 0.51, "Time = "+ "%.5f" % float(strtime))
    
    #plt.grid(True)
    plt.axes().set_aspect('equal')
    plt.tight_layout()
    #plt.savefig(PATH_FIG + "/veo_" + "%06d" % int(strstep) +".png")
    plt.show()
    plt.close()

def plot_pressure_exact(t):

    plt.figure(figsize=(5, 5))
    

    x_st = 0.0
    x_ed = math.pi
    y_st = 0.0 
    y_ed = math.pi 

    plt.xlim([x_st, x_ed])
    plt.ylim([y_st, y_ed])

    plt.xlabel(r'x')
    plt.ylabel(r'y')

    arrx = np.linspace(x_st, x_ed, 40)
    arry = np.linspace(y_st, y_ed, 40)

    X, Y = np.meshgrid(arrx, arry)
    
    P = []
    for y in arry:
        row = []
        for x in arrx:
            v = exact_u(x, y, t)
            row.append(v)
        P.append(row)

    Q = plt.contour(X, Y, P)
    plt.colorbar()

    #fu = interpolate.interp2d(arrx, arry, U, kind='linear')
    #fv = interpolate.interp2d(arrx, arry, V, kind='linear')

    #y1new = np.arange(-0.5, 0.5, 1e-2)

    #interu  = fu(0.0, y1new)
    #interv  = fv(y1new, 0.0)

    #plt.plot(interu * 0.5, y1new, color="g") 
    #plt.plot(y1new, interv * 0.5, color="g") 


    #plt.text(-0.5, 0.51, "Time = "+ "%.5f" % float(strtime))
    
    #plt.grid(True)
    plt.axes().set_aspect('equal')
    plt.tight_layout()
    #plt.savefig(PATH_FIG + "/veo_" + "%06d" % int(strstep) +".png")
    plt.show()
    plt.close()

def plot_velocity_one(strstep, strtime):
    print "Draw " + strstep + "_" + strtime  
    fnu = PATH_RESULT + "/u_" + strstep + "_" + strtime  
    fnv = PATH_RESULT + "/v_" + strstep + "_" + strtime  

    plt.figure(figsize=(5, 5))
    

    x_st = 0.0
    x_ed = math.pi
    y_st = 0.0
    y_ed = math.pi

    plt.xlim([x_st, x_ed])
    plt.ylim([y_st, y_ed])

    plt.xlabel(r'x')
    plt.ylabel(r'y')

    pd   = read.PointData(fnu)
    arrx = pd.get_coo_x()
    arry = pd.get_coo_y()

    U    = pd.get_mat_val()
    pdv  = read.PointData(fnv)
    V    = pdv.get_mat_val()

    X, Y = np.meshgrid(arrx, arry)
    M = np.hypot(U, V)

    Q = plt.quiver(X, Y, U, V, M, clim=[0.0, 1.0], cmap="bwr")
    plt.colorbar()

    fu = interpolate.interp2d(arrx, arry, U, kind='linear')
    fv = interpolate.interp2d(arrx, arry, V, kind='linear')

    y1new = np.arange(0.0, math.pi, 1e-2)

    #interu  = fu(math.pi / 2.0, y1new)
    #interv  = fv(y1new, math.pi / 2.0)

    #plt.plot(interu * 0.5, y1new, color="g") 
    #plt.plot(y1new, interv * 0.5, color="g") 

    plt.text(-0.0, math.pi + 0.01, "Time = "+ "%.5f" % float(strtime))
    
    #plt.grid(True)
    plt.axes().set_aspect('equal')
    plt.tight_layout()
    plt.savefig(PATH_FIG + "/veo_" + "%06d" % int(strstep) +".png")
    plt.close()


def plot_error(num):
    assert(num == 1 or num == 2 or num == 3)
    print "Draw error"

    plt.figure(figsize=(5, 5))
    

    x_st = 0.0
    x_ed = math.pi
    y_st = 0.0
    y_ed = math.pi

    #plt.xlim([x_st, x_ed])
    #plt.ylim([y_st, y_ed])
    plt.yscale("log")

    plt.ylabel(r'Error')
    plt.xlabel(r'Time')

    matv = read_error_file("v")
    matu = read_error_file("u") 
    matp = read_error_file("p") 

    lt = _col(matv, 1)
    lv = _col(matv, num + 1)
    lu = _col(matu, num + 1)
    lp = _col(matp, num + 1)
    lg = []
    lgt = []
    p = plt.plot(lt, lv, )
    lg.append(p[0])
    lgt.append("v")
    p = plt.plot(lt, lu, ".")
    lg.append(p[0])
    lgt.append("u")
    p = plt.plot(lt, lp)
    lg.append(p[0])
    lgt.append("p")

    plt.legend(lg, lgt)

    plt.tight_layout()
    plt.savefig(PATH_FIG + "/err"+str(num)+"_" + str(MESH) +".png")
    plt.close()

def plot_all():
    matfu = file_name(PATH_RESULT, "u")
    matfc = []
    for one in matfu:
        matfc.append(one)

    multiprocessing.freeze_support()
    pool = multiprocessing.Pool()
    cpus = multiprocessing.cpu_count() / 2
    results = []
    cmatfs = split(matfc, cpus)
    print "Thread num : ", len(cmatfs)
    for i in xrange(0, cpus):
        mat = cmatfs[i]
        for one in mat:
            result = pool.apply_async(plot_velocity_one, args=(one[1], one[2],))
            #result = pool.apply_async(plot_pressure_one, args=(one[1], one[2],))
            results.append(result)

    pool.close()
    pool.join()

    #os.system("convert -delay 20 -loop 0 "+PATH_FIG+"/veo_*.png "+PATH_FIG+"/veo.gif")
    #os.system("convert -delay 20 -loop 0 "+PATH_FIG+"/pre_*.png "+PATH_FIG+"/pre.gif")
def get_error_mat(n):
    assert(n == "u" or n == "v" or n == "p")
    matev = []
    for para in run.PARA:
        mesh = para[0]
        global MESH
        MESH = mesh
        global PATH_RESULT
        PATH_RESULT = "./res_" + str(mesh)
        mate = read_error_file(n)
        lastrow = mate[-1];
        lastrow.insert(0, para)
        matev.append(lastrow)
    return matev

def plot_all_error(n):
    assert(n == "u" or n == "v" or n == "p")
    mate = get_error_mat(n)

    print "Draw error"

    plt.figure(figsize=(5, 5))
    

    x_st = 0.0
    x_ed = math.pi
    y_st = 0.0
    y_ed = math.pi

    #plt.xlim([x_st, x_ed])
    #plt.ylim([y_st, y_ed])
    plt.yscale("log")

    plt.ylabel(r'Error')
    plt.xlabel(r'Time')

    lt = _col(mate, 0)
    le1 = _col(mate, 1)
    le2 = _col(mate, 2)
    le3 = _col(mate, 3)
    lg  = []
    lgt = []
    p = plt.plot(lt, le1, )
    lg.append(p[0])
    lgt.append("Error 1")
    p = plt.plot(lt, le2, ".")
    lg.append(p[0])
    lgt.append("Error 2")
    p = plt.plot(lt, le3)
    lg.append(p[0])
    lgt.append("Error i")

    plt.legend(lg, lgt)

    plt.tight_layout()
    plt.savefig(str(n) +".png")
    plt.close()


def main():
    for para in run.PARA:
        mesh = para[0]
        global MESH
        MESH = mesh
        global PATH_RESULT
        PATH_RESULT = "./res_" + str(mesh)
        os.system("mkdir fig_" + str(mesh))
        global PATH_FIG
        PATH_FIG = "./fig_"+str(mesh)
        #plot_all()
        plot_error(1)
        plot_error(2)
        plot_error(3)
    plot_all_error("u")
    plot_all_error("v")
    plot_all_error("p")

if __name__ == '__main__':
    main()
