#! /usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import unicode_literals, print_function
import os 
import re
import sys
import numpy as np
import math
import matplotlib
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import pylab
import argparse
import numpy as np
import string
import math
import operator
from scipy import ndimage

from PIL import Image
from PIL import ImageEnhance

import re

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = 12

class TextFile:
    def __init__(self, filename):
        self._fn = filename
        self._read_file(filename)
        self._get_dict()
        self._get_data()

    def _read_file(self, filename):
        f = open(filename, 'r')
        self._content=[]
        for i, line in enumerate(f):
            self._content.append(line)

    def is_number(self, value):
        try:
            float(value)
            return True
        except ValueError:
            return False

    def _get_dict(self):
        self._dict = {}
        for line in self._content:
            flag, string = self._parse_line_dict(line)
            if flag:
                # trim ##
                ns = string[3:-1]
                ns.strip()
                ans = ns.split(":")
                assert(len(ans) == 2)
                key = ans[0].strip()
                val = ans[1].strip()
                self._dict[key] = val

    def _get_data(self):
        self._mat = []
        for line in self._content:
            arr = self._parse_line_data(line)
            if len(arr) > 0:
                self._mat.append(arr)

    def get_config(self):
        return self._dict

    def get_data(self):
        return self._mat

    def _parse_line_data(self, line):
        # sperated by ,
        # trim space
        arr = []
        p = re.compile(r'#.+')
        if re.match(p, line):
            return arr
        
        sl = line.split(",")
        for string in sl:
            ns = string.strip()
            if self.is_number(ns):
                arr.append(float(ns))
        return arr


    def _parse_line_dict(self, line):
        # begin with ## 
        # seperated by :
        # trim space
        p = re.compile(r'##\s.+:.+')

        m = re.match(p, line)
        if m:
            return True, m.string
        else:
            return False, ""

    def get_coo_x(self):
        nx = float(self._dict["NX"])
        arr = []
        for i in range(0, int(nx)):
            arr.append(self._mat[i][0])
        return arr

    def get_coo_y(self):
        nx  = int(self._dict["NX"])
        ny  = int(self._dict["NY"])
        dim = int(self._dict["Dim"])
        arr = []
        if dim >=2:
            for i in range(0, int(ny)):
                arr.append(self._mat[i* nx][1])
            return arr
        else:
            return [0]

    def get_mat_val(self):
        assert(int(self._dict["Dim"]) ==2)
        nx  = int(self._dict["NX"])
        ny  = int(self._dict["NY"])
        dim = int(self._dict["Dim"])
        mat = []
        count = 0
        for j in range(0, ny):
            row = []
            for i in range(0, nx):
                row.append(self._mat[count][3])
                count += 1
            mat.append(row)
        return mat

    def get_arr_val(self):
        assert(int(self._dict["Dim"]) == 1)
        nx  = int(self._dict["NX"])
        dim = int(self._dict["Dim"])
        arr = []
        for i in range(0, nx):
            arr.append(self._mat[i][3])
        return arr


class PointData(TextFile):
    def __init__(self, filename):
        TextFile.__init__(self, filename)
        self.check_file()

    def check_file(self):
        assert("Size" in self._dict)


class CenterScalar:
    def __init__(self, filename):
        self._fn = filename
        self._read_file(filename)
        
    def get_filename():
        return self._fn

    def get_dim(self):
        return self._dim

    def get_nx(self):
        return self._nx

    def get_ny(self):
        return self._ny

    def get_nz(self):
        return self._nz

    def _read_file(self, filename):
        self._size = 0
        self._dim  = 0
        self._nx   = 0
        self._ny   = 0
        self._nz   = 0
        f = open(filename, 'r')
        i = 0
        self._mat_f=[]
        for i, line in enumerate(f):
            if i == 0:
                first = line.split(" ")
                assert(len(first) == 6)
                self._size = int(first[1].split(":")[1])

                self._dim = int(first[2].split(":")[1])
                self._nx = int(first[3].split(":")[1])
                self._ny = int(first[4].split(":")[1])
                self._nz = int(first[5].split(":")[1])
            else:
                lines = line.split(',')
                i=i+1
                arr_num=[];
                for numstr in lines:
                    numstr.strip()
                    arr_num.append(float(numstr))
                self._mat_f.append(arr_num)

    def get_coo_x(self):
        arrx = []
        for i in range(0, int(self._nx)):
            arrx.append(self._mat_f[i][0])
        return arrx

    def get_coo_y(self):
        arrx = []
        if self._dim >=2:
            for i in range(0, int(self._ny)):
                arrx.append(self._mat_f[i* int(self._nx)][1])
            return arrx
        else:
            return [0]

    def get_mat_val(self):
        assert(self._dim ==2)
        mat = []
        count = 0
        for j in range(0, self._ny):
            row = []
            for i in range(0, self._nx):
                row.append(self._mat_f[count][3])
                count += 1
            mat.append(row)
        return mat


def col(matrix, i):
    return [row[i] for row in matrix]

def read_residual_file(filename):
    f = open(filename, 'r')
    i = 0
    multi =[]
    _mat_f=[]
    for i, line in enumerate(f):
        line.strip()
        if line == "\n":
            multi.append(_mat_f)
            _mat_f = []
        else:
            lines = line.split(',')
            arr_num=[];
            for numstr in lines:
                numstr.strip()
                arr_num.append(float(numstr))
                _mat_f.append(arr_num)
    return multi

def test():
    print("====== start test =====")
    pd = PointData("p_0_0") 
    arr = pd.get_coo_y()
    arr = pd.get_coo_x()
    print(arr)
    print("====== end test =====")

if __name__ == '__main__':
    test()
    pass    