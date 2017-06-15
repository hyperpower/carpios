import os, sys
import copy
import numpy as np
import string
import math
import operator
import shutil

def dirs(path):
    dirs = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]
    return dirs

def clean(path, fo):
    print "clean ====== "
    files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
    dfs = []
    for f in files:
        if f not in fo:
            dfs.append(f)
    for df in dfs:
        print "Remove file -> ", df
        os.remove(df)

    dirs = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]
    for d in dirs:
        print "Remove dir -> ", d
        shutil.rmtree(d)

def main():
    path  = "./"
    pathm = os.path.abspath(path)
    ds = dirs(path)
    for d in ds:
        subd = os.path.abspath(os.path.join(path, d))
        if subd not in sys.path:
            sys.path.append(subd)
            print subd
            import run as run1
            fo = copy.copy(run1.FILE_ORIGINAL)
            clean(subd, fo)
            del run1
            sys.path.remove(subd)

        os.chdir(pathm)


if __name__ == '__main__':
    main()