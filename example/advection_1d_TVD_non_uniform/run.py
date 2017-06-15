
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
    "run.py"
]

SCHEME = [
    "QUICK",
    "VanLeer"
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
    for d in dirs:
        print "Remove dir -> ", d
        shutil.rmtree(d)

def build():
    print "cmake --------------------------------- "
    # cmake ====
    os.system("cmake .")
    print "make  --------------------------------- "
    os.system("make")
    print "run   --------------------------------- "
    for s in SCHEME:
        os.system("./build/main " + s)
    print "plot   -------------------------------- "
    os.system("python plot.py")

def main():
    clean() 
    os.system("mkdir result")
    os.system("mkdir fig")
    build()

if __name__ == '__main__':
    main()
