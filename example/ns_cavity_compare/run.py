
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
    "gerris_result",
    "monitor.py",
    "save"
]

PARA = [
    [1000, 32, "upwind2"]
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
        re   = para[0]
        mesh = para[1]
        adv  = para[2]
        os.system("mkdir res_"+str(re)+"_"+str(mesh)+"_"+str(adv))

    for para in PARA:
        re   = para[0]
        mesh = para[1]
        adv  = para[2]
        os.system("./build/main "+str(re)+" "+str(mesh)+" "+str(adv))
    print "plot   -------------------------------- "
    #os.system("python plot.py")

def main():
    clean() 
    os.system("mkdir result")
    os.system("mkdir fig")
    build()

if __name__ == '__main__':
    main()
