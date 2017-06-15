
import os
import numpy as np
import string
import math
import operator
import shutil


FILE_ORIGINAL = [
    "main.cpp",
    "CMakeLists.txt",
    "plot.py",
    "run.py",
    "monitor.py",
    "save"
]

PARA = [
    [100],
    [200],
    [400],
    [800],
    [1600],
]

def clean():
    print "clean ====== "
    path = "./"
    files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
    dfs = []
    for f in files:
        if f not in FILE_ORIGINAL:
            dfs.append(f)
    for df in dfs:
        print "Remove file -> ", df
        os.remove(df)

    dirs = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]
    dfs  = []
    for f in dirs:
        if f not in FILE_ORIGINAL:
            dfs.append(f)
    for d in dfs:
        print "Remove dir -> ", d
        shutil.rmtree(d)

def build():
    print "cmake --------------------------------- "
    # cmake ====
    os.system("cmake .")
    print "make  --------------------------------- "
    os.system("make")
    print "run   --------------------------------- "
    for para in PARA:
        mesh   = para[0]
        os.system("mkdir res_"+str(mesh))

    for para in PARA:
        mesh   = para[0]
        os.system("./build/main "+str(mesh))
    print "plot   -------------------------------- "
    #os.system("python plot.py")

def main():
    clean() 
    build()

if __name__ == '__main__':
    main()
