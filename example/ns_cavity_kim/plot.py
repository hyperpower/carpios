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

PATH_EXAMPLE  = os.path.abspath(os.path.join(__file__, "../"))
PATH_THIS     = os.path.abspath(__file__)
PATH_RESULT   = os.path.abspath("./result")
PATH_FIG      = os.path.abspath("./fig")
PATH_PROJECT  = os.path.abspath(os.path.join(PATH_EXAMPLE, "../.."))
PATH_PYSCRIPT = os.path.abspath(os.path.join(PATH_PROJECT, "pyscript"))

sys.path.append(PATH_PYSCRIPT)
import read

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

def get_gfs_result_u(path, re, mesh):
    filename  = path + "/xprof_" + str(re) + "_" + str(mesh)

    # If you need to open a file instead:
    res = []
    f = open(filename)
    for line in f:
        fields = line.strip().split()
        # Array indices start at 0 unlike AWK
        if fields[2] != '2:x':
            res.append([float(fields[2]), float(fields[6])])
    return res 

def get_gfs_result_v(path, re, mesh):
    filename  = path + "/yprof_" + str(re) + "_" + str(mesh)

    # If you need to open a file instead:
    res = []
    f = open(filename)
    for line in f:
        fields = line.strip().split()
        # Array indices start at 0 unlike AWK
        if fields[1] != '1:t':
            res.append([float(fields[1]), float(fields[7])])
    return res 

def get_gerris_result(path, re, level, uorv):
    if uorv == "u":
        return get_gfs_result_u(path, re, level)
    elif uorv == "v":
        return get_gfs_result_v(path, re, level)



def plot_velocity_one(strstep, strtime):
    print "Veo" , strstep," ", strtime
    fnu = PATH_RESULT + "/u_" + strstep + "_" + strtime  
    fnv = PATH_RESULT + "/v_" + strstep + "_" + strtime  

    plt.figure(figsize=(5, 5))
    

    x_st = -0.5
    x_ed = 0.5
    y_st = -0.5
    y_ed = 0.5

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

    Q = plt.quiver(X, Y, U, V, M, clim=[0,1], cmap="bwr")
    plt.colorbar()

    fu = interpolate.interp2d(arrx, arry, U, kind='linear')
    fv = interpolate.interp2d(arrx, arry, V, kind='linear')

    y1new = np.arange(-0.5, 0.5, 1e-2)

    interu  = fu(0.0, y1new)
    interv  = fv(y1new, 0.0)

    plt.plot(interu * 0.5, y1new, color="g") 
    plt.plot(y1new, interv * 0.5, color="g") 

    # gerris compare
    gmat = get_gerris_result("./gerris_result", RE, 6, "u")
    gx   = np.array(_col(gmat, 0))
    gv   = np.array(_col(gmat, 1))
    plt.plot(gv * 0.5, gx, color="k") 
    gmat = get_gerris_result("./gerris_result", RE, 6, "v")
    gx   = np.array(_col(gmat, 0))
    gv   = np.array(_col(gmat, 1))
    plt.plot(gx, gv *  0.5, color="k") 
        #centeru = plt.plot(X, interu[:][0], cmap="bwr")

    plt.text(-0.5, 0.51, "Time = "+ "%.5f" % float(strtime))
    
    #plt.grid(True)
    plt.axes().set_aspect('equal')
    plt.tight_layout()
    plt.savefig(PATH_FIG + "/veo_" + "%06d" % int(strstep) +".png")
    plt.close()

def plot_pressure_one(strstep, strtime):
    print "Pre" , strstep," ", strtime
    fnp = PATH_RESULT + "/p_" + strstep + "_" + strtime  

    plt.figure(figsize=(5, 5))
    
    x_st = -0.5
    x_ed = 0.5
    y_st = -0.5
    y_ed = 0.5

    plt.xlim([x_st, x_ed])
    plt.ylim([y_st, y_ed])

    plt.xlabel(r'x')
    plt.ylabel(r'y')

    pd   = read.PointData(fnp)
    arrx = pd.get_coo_x()
    arry = pd.get_coo_y()

    P    = pd.get_mat_val()

    X, Y = np.meshgrid(arrx, arry)

    Q = plt.contour(X, Y, P, clim=[-0.1, 1.0], cmap="bwr")
    plt.colorbar()

    fp = interpolate.interp2d(arrx, arry, P, kind='linear')

    y1new = np.arange(-0.5, 0.5, 1e-2)

    interp1  = fp(0.0, y1new)
    interp2  = fp(y1new, 0.0)

    plt.plot(interp1 * 0.5, y1new, color="g") 
    plt.plot(y1new, interp2 * 0.5, color="g") 
        #centeru = plt.plot(X, interu[:][0], cmap="bwr")

    plt.text(-0.5, 0.51, "Time = "+ "%.5f" % float(strtime))
    
    #plt.grid(True)
    plt.axes().set_aspect('equal')
    plt.tight_layout()
    plt.savefig(PATH_FIG + "/pre_" + "%06d" % int(strstep) +".png")
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
            result = pool.apply_async(plot_pressure_one, args=(one[1], one[2],))
            results.append(result)

    pool.close()
    pool.join()

    #os.system("convert -delay 20 -loop 0 "+PATH_FIG+"/veo_*.png "+PATH_FIG+"/veo.gif")
    #os.system("convert -delay 20 -loop 0 "+PATH_FIG+"/pre_*.png "+PATH_FIG+"/pre.gif")

def main():
    for para in run.PARA:
        re   = para[0]
        mesh = para[1]
        adv  = para[2]
        global RE
        global MESH
        global ADV
        RE   = re
        MESH = mesh
        ADV  = adv
        global PATH_RESULT
        PATH_RESULT = "./res_"+str(re)+"_"+str(mesh)+"_"+str(adv)
        os.system("mkdir fig_"+str(re)+"_"+str(mesh)+"_"+str(adv))
        global PATH_FIG
        PATH_FIG = "./fig_"+str(re)+"_"+str(mesh)+"_"+str(adv)
        plot_all()
    #plot_velocity_one("54","0.00265886") 

if __name__ == '__main__':
    main()
