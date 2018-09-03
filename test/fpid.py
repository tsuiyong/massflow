# -*- coding: UTF-8 -*-
import h5py
import numpy as np
import pandas as pd


a = np.array([])
b = np.array([])
f = open("demo.txt")

while 1:
    lines = f.readlines(1247)
    if not lines:
        break
    for line in lines:
        print(line)
        a = np.append(a, int(line))
f.close()

f = h5py.File("FCDB.hdf5", "a")

d = f["/NuclInfo/NuclID"]

for i in range(0, 1693):
    b = np.append(b, d[i])

c = np.intersect1d(a, b)     #裂变产物核素

ds = f.create_dataset("Fission_Product_ID", (len(c),), dtype = np.int32)

for i in range(0, len(c)):
    ds[i] = c[i].astype(np.int32)

f.close()