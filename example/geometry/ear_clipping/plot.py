import matplotlib
import matplotlib.pyplot as plt
import os, sys
import numpy as np
import string
import math
import operator


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
PATH_PROJECT  = os.path.abspath(os.path.join(PATH_EXAMPLE, "../../.."))
PATH_PYSCRIPT = os.path.abspath(os.path.join(PATH_PROJECT, "pyscript"))

sys.path.append(PATH_PYSCRIPT)
import read

import run

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




def plot_one(fn, idx):
    print "Draw : ", fn, " ", idx

    plt.figure(figsize=(6, 4))

    """
    Set labels
    """
    plt.xlabel(r'x')
    plt.ylabel(r'y')

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
    filename = "./res_" + fn + "/" + "%03d" % idx 

    f = read.TextFile(filename)
    d    = f.get_config()
    data = f.get_data()

    data.append(data[0])

    x = _col(data, 0);
    y = _col(data, 1);

    plt.plot(x, y, "r", linewidth = 2)

    # draw other
    for ip in range(0, idx):
        filename = "./res_" + fn + "/" + "%03d" % ip 
        f        = read.TextFile(filename)
        data     = f.get_data()
        data.append(data[0])
        x = _col(data, 0);
        y = _col(data, 1);
        plt.plot(x, y, "b", linewidth = 1)

    # plot origin
    filename = "./" + fn 
    f        = read.TextFile(filename)
    data     = f.get_data()
    data.append(data[0])
    x = _col(data, 0);
    y = _col(data, 1);
    plt.plot(x, y, "k", linewidth = 1)


    #plt.legend(llg, scheme, loc= 'upper right')

    #plt.grid(Tru
    #plt.axes().set_aspect('equal')
    plt.tight_layout()
    plt.savefig("./res_" + fn + "/" + "%03d" % idx  +".png")
    plt.close()
    #plt.show()




def plot_all():
    matfu = file_name(PATH_RESULT, "phi")
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
            result = pool.apply_async(plot_one, args=(one[2], one[3],))
            results.append(result)

    pool.close()
    pool.join()

    os.system("convert -delay 5 -loop 0 ./fig/comp_*.png comp.gif")


def main(fn):
    # fn      = "man"
    namedir = "./res_" + fn + "/" 
    res = []
    files = [f for f in os.listdir(namedir) if os.path.isfile(os.path.join(namedir, f))]
    for f in files:
        if len(f) == 3:
            res.append(int(f))

    res.sort()

    for idx in res:
        plot_one(fn, idx)

    os.system("convert -delay 10 -loop 0 " + namedir + "/*.png "+fn+".gif")

if __name__ == '__main__':
    for n in run.DATA_FILES:
        main(n)


