#!/usr/bin/python
# -*- coding:utf-8 -*-

import pandas as pd

df = pd.read_table('input.txt', sep='\s+')

for i in range(4):
    print(df.iloc[i, 1])

