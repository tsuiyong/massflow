#!/usr/bin/python
# -*- coding:utf-8 -*-

# 衰变阶段计算

    n0 = np.zeros(1693, dtype=np.float64)  # 置零
    n0[nuclid] = mass3[year]

    # 计算初始时刻的值
    y, y_act, y_q, y_toxi = CRAM48(n0, 0)
    concentration = np.array([])
    concentration = np.append(concentration, y[nuclid])
    activity = np.array([])
    activity = np.append(activity, y_act[nuclid])
    heat = np.array([])
    heat = np.append(heat, y_q[nuclid])
    toxicity = np.array([])
    toxicity = np.append(toxicity, y_toxi[nuclid])

    tl = input()
    for i in range(1, int(tl) + 1, int(int(tl) * 0.01)):
        dt = float(i) * 365.25 * 24. * 3600.
        y, y_act, y_q, y_toxi = CRAM48(n0, dt)
        concentration = np.append(concentration, y[nuclid])
        activity = np.append(activity, y_act[nuclid])
        heat = np.append(heat, y_q[nuclid])
        toxicity = np.append(toxicity, y_toxi[nuclid])

    f.close()