#!/usr/bin/env python3
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
    if len(s) == 5:
        kg = mole * int(s[2:4]) / 1000
    else:
        if s[2] == 0:
            kg = mole * int(s[3:5]) / 1000
        else:
            kg = mole * int(s[2:5]) / 1000
    f.close()
    return kg


# 卸堆裂变产物质量计算
def fission_product(power):
    year = 109
    Mass_fp = 0
    min, max, PowerDensity = initialize(power)
    fp = [131, 132, 140, 141, 142, 143, 152, 153, 154, 155, 156, 157, 158, 170, 171, 172, 173, 174, 175, 176,
          177, 178, 179, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 207, 208, 209, 210,
          211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 228, 229, 230, 231, 232, 233,
          234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253,
          254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273,
          274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 293, 294,
          295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 308, 309, 310, 311, 312, 313, 314, 315,
          316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335,
          336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355,
          356, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376,
          377, 378, 379, 380, 381, 382, 383, 384, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397,
          398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 417, 418,
          419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438,
          440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459,
          460, 461, 462, 463, 464, 465, 467, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481,
          482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 498, 499, 500, 502, 503,
          504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523,
          524, 525, 527, 529, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546,
          547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566,
          567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 586,
          587, 588, 589, 590, 592, 593, 594, 595, 596, 597, 598, 599, 600, 601, 602, 603, 604, 605, 606, 607,
          608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627,
          628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 647, 648, 649,
          650, 651, 652, 653, 654, 655, 656, 657, 658, 659, 660, 661, 662, 663, 664, 665, 666, 667, 668, 669,
          670, 671, 672, 673, 674, 675, 676, 677, 678, 679, 681, 682, 683, 684, 685, 686, 687, 688, 689, 690,
          691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710,
          711, 712, 713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 723, 724, 725, 726, 727, 728, 729, 730,
          731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 745, 746, 747, 748, 749, 750,
          751, 752, 753, 754, 755, 756, 757, 758, 759, 760, 761, 762, 763, 764, 765, 766, 767, 768, 769, 770,
          771, 772, 773, 774, 775, 777, 778, 779, 780, 781, 782, 783, 784, 785, 786, 787, 788, 790, 791, 792,
          793, 794, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 805, 806, 807, 808, 809, 810, 811, 812,
          813, 814, 815, 816, 817, 818, 819, 820, 821, 822, 823, 824, 825, 826, 827, 828, 829, 830, 831, 832,
          833, 834, 835, 836, 837, 838, 839, 840, 841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852,
          853, 854, 855, 856, 857, 858, 859, 860, 861, 862, 863, 864, 865, 866, 867, 868, 869, 870, 871, 872,
          873, 874, 875, 876, 877, 878, 879, 880, 881, 882, 883, 884, 885, 886, 887, 888, 889, 890, 891, 892,
          894, 895, 896, 897, 898, 899, 900, 901, 902, 903, 904, 905, 906, 907, 908, 909, 910, 912, 913, 914,
          915, 916, 917, 918, 919, 920, 921, 922, 923, 924, 925, 926, 927, 928, 929, 930, 931, 932, 933, 934,
          935, 936, 937, 938, 939, 940, 941, 942, 944, 946, 947, 948, 949, 950, 951, 952, 953, 954, 955, 956,
          957, 958, 959, 960, 961, 962, 963, 964, 965, 966, 967, 968, 969, 970, 971, 974, 975, 976, 977, 978,
          979, 980, 981, 982, 983, 984, 985, 986, 987, 988, 989, 990, 991, 992, 993, 994, 995, 996, 997, 998,
          999, 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1008, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017,
          1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034,
          1036, 1038, 1040, 1041, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1053, 1054, 1055,
          1056, 1057, 1058, 1059, 1060, 1061, 1062, 1063, 1064, 1065, 1066, 1068, 1070, 1071, 1072, 1073, 1074,
          1075, 1076, 1077, 1078, 1079, 1080, 1081, 1082, 1083, 1084, 1085, 1086, 1087, 1088, 1089, 1091, 1093,
          1095, 1096, 1097, 1098, 1099, 1100, 1101, 1102, 1103, 1104, 1105, 1106, 1107, 1108, 1109, 1110, 1111,
          1112, 1113, 1114, 1115, 1116, 1117, 1118, 1119, 1120, 1121, 1122, 1123, 1124, 1125, 1126, 1127, 1128,
          1129, 1130, 1131, 1132, 1133, 1134, 1135, 1136, 1137, 1138, 1139, 1140, 1141, 1143, 1144, 1145, 1146,
          1147, 1148, 1149, 1150, 1151, 1152, 1153, 1154, 1155, 1156, 1157, 1158, 1159, 1160, 1161, 1162, 1163,
          1164, 1165, 1166, 1167, 1168, 1169, 1170, 1171, 1172, 1173, 1174, 1175, 1176, 1177, 1178, 1179, 1180,
          1181, 1182, 1183, 1184, 1185, 1186, 1187, 1188, 1189, 1190, 1191, 1192, 1193, 1194, 1195, 1196, 1197,
          1198, 1199, 1200, 1201, 1202, 1203, 1204, 1205, 1206, 1207, 1208, 1209, 1210, 1211, 1214, 1216, 1219,
          1220, 1221, 1222, 1223, 1224, 1225, 1226, 1227, 1228, 1229, 1230, 1231, 1232, 1233, 1234, 1235, 1236,
          1237, 1239, 1241, 1243, 1244, 1245, 1246, 1247, 1249, 1250, 1251, 1252, 1253, 1254, 1255, 1256, 1257,
          1258, 1259, 1260, 1261, 1262, 1263, 1265, 1268, 1269, 1270, 1271, 1272, 1273, 1274, 1275, 1276, 1277,
          1278, 1279, 1281, 1282, 1283, 1284, 1285, 1286, 1287, 1288, 1290, 1291, 1292, 1293, 1294, 1295, 1296,
          1297, 1298, 1299, 1300, 1301, 1302, 1303, 1304, 1305, 1306, 1307, 1308, 1309, 1310, 1313, 1314, 1315,
          1316, 1317, 1318, 1319, 1320, 1321, 1322, 1323, 1324, 1325, 1326, 1327, 1328, 1329, 1330, 1332, 1333,
          1334, 1335, 1336, 1337, 1338, 1339, 1340, 1341, 1342, 1343, 1344, 1345, 1346, 1347, 1348, 1349, 1350,
          1352, 1353, 1355, 1356, 1357, 1358, 1359, 1360, 1361, 1362, 1369, 1370, 1372, 1373, 1374, 1375, 1385,
          1386]
    for j in fp:
        m1 = mass_off(year, min, j)
        m2 = mass_off(year, max, j)
        x = [PowerDensity[min], PowerDensity[max]]
        y = [m1, m2]
        func = lagrange(x, y)
        m3 = func(float(power))
        Mass_fp = Mass_fp + mole2kg(m3, j)
    return Mass_fp


# 卸堆超铀元素质量计算
def transuranium(power):
    year = 109
    Mass_tru = 0
    min, max, PowerDensity = initialize(power)
    tru = [1623, 1624, 1625, 1626, 1627, 1628, 1629, 1630, 1631, 1632, 1633, 1634, 1635, 1636, 1637, 1638, 1639,
           1640, 1641, 1642, 1643, 1644, 1645, 1646, 1647, 1648, 1649, 1650, 1651, 1652, 1653, 1654, 1655, 1656,
           1657, 1658, 1659, 1660, 1661, 1662, 1663, 1664, 1665, 1666, 1667, 1668, 1669]
    for i in tru:
        m1 = mass_off(year, min, i)
        m2 = mass_off(year, max, i)
        x = [PowerDensity[min], PowerDensity[max]]
        y = [m1, m2]
        func = lagrange(x, y)
        m3 = func(float(power))
        Mass_tru = Mass_tru + mole2kg(m3, i)
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
    return mass3[0]  # 初装量


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

    if decay_time == 100:
        time = [t for t in range(1, decay_time + 1, int(decay_time / 100))]
    else:
        time = [t for t in range(1, decay_time + 1, int(decay_time / 1000))]
    act = [radio_activity(t, mass_off) for t in time]
    plt.loglog(time, act)
    plt.savefig(sys.path[0] + "/figs/" + "activity_" + str(decay_time) + ".png")
    plt.show()
    return


def required_nuclides(mass_off):
    f = h5py.File(sys.path[0] + "/lib/FCDB.hdf5", "r")
    namedset = f["/NuclInfo/NuclName"]
    file = open(sys.path[0] + "/nuclide.txt", "w")
    nuclides = [1598, 1612, 1613, 1614, 1615, 1616, 1617, 1618, 1619, 1627, 1636, 1637, 1638, 1639,
                1640, 1648, 1650, 1651, 1660, 1661, 1662, 1663, 1664]
    for i in nuclides:
        if i == 1615:
            file.write(namedset[i].decode('utf-8') + " " + str((mole2kg(mass_off[i], i)
                                                                + mole2kg(mass_off[i + 1], i + 1)) * 1000) + "\n")
        elif i == 1616:
            continue
        else:
            file.write(namedset[i].decode('utf-8') + " " + str(mole2kg(mass_off[i], i) * 1000) + "\n")
    file.close()
    f.close()
