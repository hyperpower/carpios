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

def main():
    path  = "./"
    pathm = os.path.abspath(path)
    ds = dirs(path)
    for d in ds:
        subd = os.path.abspath(os.path.join(path, d))
        os.chdir(subd)
        os.system("python run.py")
        os.chdir(pathm)


if __name__ == '__main__':
    main()