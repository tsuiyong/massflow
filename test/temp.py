#!/usr/bin/python
# -*- coding:utf-8 -*-

from func import hb


decay_time = 1000000
s1 = [t for t in range(1, 1001)]
s2 = [t for t in range(1001, decay_time + 1, int(decay_time / 1000))]

time = s1 + s2
print(time)