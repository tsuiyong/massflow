#!/usr/bin/python
# -*- coding:utf-8 -*-

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib  # 用于全局参数设置
import os
import sys
from scipy.interpolate import lagrange  # Lagrange 插值方法
from decay.decay import CRAM48  # 导入decay模块


# 卸堆核素浓度计算
def mass_off(time, index1, index2):
    f = h5py.File(sys.path[0] + "/lib/FCDB.hdf5", "r")
    s1 = "MSBR_0" + str(index1 + 1) + "/"
    s2 = str(time) + "/CoreInventory"
    dset = f[os.path.join("/MSR/", s1, "TimeEvolution/", s2)]
    concentration = dset[index2]
    f.close()
    return concentration


# 程序初始化，对功率密度进行插值
def initialize(P):
    PowerDensity = np.array([])
    num = 3  # 功率密度个数
    f = h5py.File(sys.path[0] + "/lib/FCDB.hdf5", "r")
    for i in range(1, num + 1):
        s1 = "MSBR_0" + str(i)
        grp = f[os.path.join("/MSR/", s1)]
        for key, value in grp.attrs.items():
            if key == "SpecificPower":
                PowerDensity = np.append(PowerDensity, value)

    # 根据给定的功率密度对核素浓度进行插值
    for i in range(num):
        if float(P) == PowerDensity[i]:
            pass  # 有待添加
    for i in range(num - 1):
        if PowerDensity[i] < float(P) < PowerDensity[i + 1]:
            min0 = i  # mass 函数输入参数
            max0 = i + 1  # mass 函数输入参数

    if float(P) < PowerDensity[0] or float(P) > PowerDensity[num - 1]:
        print("Out of the boundary!!")
        exit()
    f.close()
    return min0, max0, PowerDensity


# 根据核素符号查询核素序号
def alias2id(alias):
    f = h5py.File(sys.path[0] + "/lib/FCDB.hdf5", "r")
    dset = f["/NuclInfo/NuclName"]  # 打开元素符号数据集
    for i in range(0, 1693):
        alias0 = dset[i].decode('utf-8')  # 将bytes数据转换为str类型
        if alias == alias0:
            nuclid = i  # 核素序号
            break
    f.close()
    return nuclid


# 物质的量(mole)转换为质量(Kg)
def mole2kg(mole, id):
    f = h5py.File(sys.path[0] + "/lib/FCDB.hdf5", "r")
    dataset = f["/NuclInfo/NuclID"]
    s = str(dataset[id])
    kg = mole * int(s[2:5]) / 1000
    f.close()
    return kg


# 卸堆裂变产物质量计算
def fission_product(power):
    year = 109
    Mass_fp = 0
    min, max, PowerDensity = initialize(power)
    f = h5py.File(sys.path[0] + "/lib/FCDB.hdf5", "r")
    dataset = f["Fission_Product_ID"]
    idset = f["/NuclInfo/NuclID"]
    for i in range(0, len(dataset)):
        for j in range(0, 1693):
            if dataset[i] == idset[j]:
                m1 = mass_off(year, min, j)
                m2 = mass_off(year, max, j)
                x = [PowerDensity[min], PowerDensity[max]]
                y = [m1, m2]
                func = lagrange(x, y)
                m3 = func(float(power))
                break
        s = str(dataset[i])
        print(m3)
        Mass_fp = Mass_fp + m3 * int(s[2:5])

    Mass_fp = Mass_fp / 1000
    f.close()
    return Mass_fp


# 卸堆超铀元素质量计算
def transuranium(power):
    year = 109
    Mass_tru = 0
    min, max, PowerDensity = initialize(power)
    f = h5py.File(sys.path[0] + "/lib/FCDB.hdf5", "r")
    dataset = f["TRU_ID"]
    idset = f["/NuclInfo/NuclID"]
    for i in range(0, len(dataset)):
        for j in range(0, 1693):
            if dataset[i] == idset[j]:
                m1 = mass_off(year, min, j)
                m2 = mass_off(year, max, j)
                x = [PowerDensity[min], PowerDensity[max]]
                y = [m1, m2]
                func = lagrange(x, y)
                m3 = func(float(power))
                break
        s = str(dataset[i])
        print(m3)
        Mass_tru = Mass_tru + m3 * int(s[2:5])

    Mass_tru = Mass_tru / 1000
    f.close()
    return Mass_tru


# 核素随时间演化计算
def mass(time, index1, index2):
    f = h5py.File(sys.path[0] + "/lib/FCDB.hdf5", "r")
    mass_flow = np.array([])
    for i in range(0, time + 1):
        s1 = "MSBR_0" + str(index1 + 1) + "/"
        s2 = str(i) + "/CoreInventory"
        dset = f[os.path.join("/MSR/", s1, "TimeEvolution/", s2)]
        mass_flow = np.append(mass_flow, dset[index2])
    f.close()
    return mass_flow


# 单一核素初装浓度计算
def nuc1(power, nuclid):
    year = 109  # 年份，单位：年
    min, max, PowerDensity = initialize(power)
    # 准备插值数据
    mass1 = mass(year, min, nuclid)
    mass2 = mass(year, max, nuclid)
    # 插值
    mass3 = np.array([])
    for i in range(0, year + 1):
        x = [PowerDensity[min], PowerDensity[max]]
        y = [mass1[i], mass2[i]]
        func = lagrange(x, y)
        mass3 = np.append(mass3, func(float(power)))
    return mass3[0]# 初装量


# 单一核素卸堆浓度计算
def nuc2(power, nuclid):
    year = 109  # 年份，单位：年
    min, max, PowerDensity = initialize(power)
    # 准备插值数据
    mass1 = mass(year, min, nuclid)
    mass2 = mass(year, max, nuclid)
    # 插值
    mass3 = np.array([])
    for i in range(0, year + 1):
        x = [PowerDensity[min], PowerDensity[max]]
        y = [mass1[i], mass2[i]]
        func = lagrange(x, y)
        mass3 = np.append(mass3, func(float(power)))
    return mass3[109]  # 卸堆量


# 所有核素衰变计算
def decay_cal(t, mass_off):
    n0 = mass_off
    dt = t * 365.25 * 24. * 3600.
    y, y_act, y_q, y_toxi = CRAM48(n0, dt)
    return y, y_act, y_q, y_toxi


# 反应堆运行阶段核素浓度演化绘图
def paint_core(time, Alias, power, massflow):
    matplotlib.rcParams['font.family'] = 'sans-serif'  # 绘图风格
    matplotlib.rcParams['font.sans-serif'] = 'Times New Roman'  # 字体设置Times New Roman

    t = np.arange(0, time + 1)
    plt.title('Mass flow of ' + Alias + "(Power Density =" + str(power) + ")", \
              fontsize=20)
    plt.xlabel('Time(year)', fontsize=15)
    plt.ylabel('Concentration of nuclide($m^{-3}$)', fontsize=15)
    plt.grid(True)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.plot(t, massflow)
    plt.savefig(sys.path[0] + "/figs/" + "Massflow_" + Alias + ".png")
    plt.show()
    return


# 衰变阶段参数绘图
def paint_decay(time, name, value, kind):
    matplotlib.rcParams['font.family'] = 'sans-serif'  # 绘图风格
    matplotlib.rcParams['font.sans-serif'] = 'Times New Roman'  # 字体设置Times New Roman

    if kind == 'Concentration':
        unit = 'Concentration(mol)'
    elif kind == 'Activity':
        unit = 'Activity(Ci)'
    elif kind == 'Decay heat':
        unit = 'Decay heat(W)'
    else:
        unit = 'Toxicity(Sv)'

    t = np.arange(1, time + 1, int(time * 0.01))
    plt.title(kind + ' of ' + name, fontsize=20)
    plt.xlabel('Time(year)', fontsize=15)
    plt.ylabel(unit, fontsize=15)
    plt.grid(True)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.plot(t, value)
    plt.savefig(sys.path[0] + "/figs/" + kind + "_" + name + ".png")
    plt.show()
    return


def nuclide_evolution(power, nuclid):
    year = 109  # 年份，单位：年
    min, max, PowerDensity = initialize(power)
    # 准备插值数据
    mass1 = mass(year, min, nuclid)
    mass2 = mass(year, max, nuclid)
    # 插值
    mass3 = np.array([])
    for i in range(0, year + 1):
        x = [PowerDensity[min], PowerDensity[max]]
        y = [mass1[i], mass2[i]]
        func = lagrange(x, y)
        mass3 = np.append(mass3, func(float(power)))
    return mass3


def radio_activity(time, mass_off):
    concentration, activity, heat, toxicity = decay_cal(time, mass_off)
    return np.sum(activity)


def paint_activity(decay_time, mass_off):
    matplotlib.rcParams['font.family'] = 'sans-serif'  # 绘图风格
    matplotlib.rcParams['font.sans-serif'] = 'Times New Roman'  # 字体设置Times New Roman

    plt.title("Radio activity" + ' of ' + str(decay_time) + "years", fontsize=20)
    plt.xlabel('Time(year)', fontsize=15)
    plt.ylabel("Total activity(Ci)", fontsize=15)
    plt.grid(True)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    time = [t for t in range(1, decay_time + 1, int(decay_time / 100))]
    act = [radio_activity(t, mass_off) for t in time]
    plt.loglog(time, act)
    plt.savefig(sys.path[0] + "/figs/" + "activity_" + str(decay_time) + ".png")
    plt.show()
    return


