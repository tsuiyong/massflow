# encoding: UTF-8
import numpy as np
from decay import CRAM48 

# 衰变算例 1mol Np-237 衰变100万年
n0=np.zeros(1693,dtype=np.float64)
n0[1627]=1.0 #9.6963e-1
dt= 1000000.0 * 365.25 * 24. * 3600.

# 数组分别为时刻末浓度、放射性活度、衰变热以及放射性毒性
y, y_act, y_q, y_toxi = CRAM48(n0, dt)

print(y[1627])
print(y_act[1627])
print(y_q[1627])
print(y_toxi[1627])