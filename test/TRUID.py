# -*- coding: UTF-8 -*-

import h5py
import numpy as np

f = h5py.File("FCDB.hdf5", "a")
d = f["/NuclInfo/NuclID"]

TRU = np.array([])
for i in range(0, len(d)):
    s = str(d[i])
    if len(s) > 5:
        if int(s[0:2]) >= 93 and int(s[0:2]) <= 96:
            TRU = np.append(TRU, d[i])

ds = f.create_dataset("TRU_ID", (len(TRU),), dtype=np.int32)

for i in range(0, len(TRU)):
    ds[i] = TRU[i].astype(np.int32)
