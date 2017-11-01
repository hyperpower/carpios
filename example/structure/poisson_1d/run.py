from __future__ import print_function 
import time
import platform
import os
import numpy as np
import string
import math
import operator
import shutil
import cpuinfo


FILE_ORIGINAL = [
    "main.cpp",
    "CMakeLists.txt",
    "plot.py",
    "run.py",
    "plotly_draw.py"
]

def clean():
    print("clean ====== ")
    path = "./"
    files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
    dfs = []
    for f in files:
        if f not in FILE_ORIGINAL:
            dfs.append(f)
    for df in dfs:
        print("Remove file -> %s" % df)
        os.remove(df)

    dirs = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]
    for d in dirs:
        print("Remove dir -> %s" % d)
        shutil.rmtree(d)

def build():
    print("cmake --------------------------------- ")
    # cmake ====
    os.system("cmake .")
    print("make  --------------------------------- ")
    os.system("make")
    print("run   --------------------------------- ")
    ci = cpuinfo.get_cpu_info()
    print("System     : %s" % platform.system())
    print("CPU        : %s" % ci["brand"])
    t0 = time.time()
    os.system("./build/main")
    t1 = time.time()
    total = t1 - t0
    print("Time spend : %10.5e s" % total)
    print("plot   -------------------------------- ")
    os.system("python3 plot.py")
    os.system("python3 plotly_draw.py")


def main():
    clean() 
    build()

if __name__ == '__main__':
    main()
