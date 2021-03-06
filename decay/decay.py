# encoding: UTF-8

"""
Author: ShaoPeng Xia @ DRP.SINAP
Date: 11:37 2018-6-26

Description
===========

This script is used to calculate the mass flow for fuel cycle code TMSR-NUDEAS.
"""

import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

def CRAM48(n0, dt):
    """ 48阶切比雪夫有理近似求解衰变方程，且采用数值更为稳定的IPF格式

    输入参数
    ----------
    n0 : np.array
        初始核素浓度，单位 mol
    dt : float
        衰变时间，单位 s

    返回结果
    -------
    y : np.array
        衰变时刻末的核素浓度，单位 mol
    y_act : np.array
        衰变时刻末的核素放射性活度，单位 Ci
    y_q : np.array
        衰变时刻末的衰变热，单位 Watts
    y_toxi : np.array
        衰变时刻末的放射性毒性，单位 Sv

    """

    #A : sp.csr_matrix 衰变系数矩阵.
    A = sp.load_npz(sys.path[0] + '/lib/Matrix_A.npz')

    #Coeff : np.array 存储各核素的衰变相关系数，包括衰变常数("Lambda")，衰变热系数("DecayHeat")以及放射性毒性系数("Toxicity")
    Coeff = np.load(sys.path[0] + '/lib/DecayCoeff.npz')

    theta_r = np.array([-4.465731934165702e+1, -5.284616241568964e+0,
                        -8.867715667624458e+0, +3.493013124279215e+0,
                        +1.564102508858634e+1, +1.742097597385893e+1,
                        -2.834466755180654e+1, +1.661569367939544e+1,
                        +8.011836167974721e+0, -2.056267541998229e+0,
                        +1.449208170441839e+1, +1.853807176907916e+1,
                        +9.932562704505182e+0, -2.244223871767187e+1,
                        +8.590014121680897e-1, -1.286192925744479e+1,
                        +1.164596909542055e+1, +1.806076684783089e+1,
                        +5.870672154659249e+0, -3.542938819659747e+1,
                        +1.901323489060250e+1, +1.885508331552577e+1,
                        -1.734689708174982e+1, +1.316284237125190e+1])
    theta_i = np.array([+6.233225190695437e+1, +4.057499381311059e+1,
                        +4.325515754166724e+1, +3.281615453173585e+1,
                        +1.558061616372237e+1, +1.076629305714420e+1,
                        +5.492841024648724e+1, +1.316994930024688e+1,
                        +2.780232111309410e+1, +3.794824788914354e+1,
                        +1.799988210051809e+1, +5.974332563100539e+0,
                        +2.532823409972962e+1, +5.179633600312162e+1,
                        +3.536456194294350e+1, +4.600304902833652e+1,
                        +2.287153304140217e+1, +8.368200580099821e+0,
                        +3.029700159040121e+1, +5.834381701800013e+1,
                        +1.194282058271408e+0, +3.583428564427879e+0,
                        +4.883941101108207e+1, +2.042951874827759e+1])
    theta = np.array(theta_r + theta_i * 1j, dtype=np.complex128)

    alpha_r = np.array([+6.387380733878774e+2, +1.909896179065730e+2,
                        +4.236195226571914e+2, +4.645770595258726e+2,
                        +7.765163276752433e+2, +1.907115136768522e+3,
                        +2.909892685603256e+3, +1.944772206620450e+2,
                        +1.382799786972332e+5, +5.628442079602433e+3,
                        +2.151681283794220e+2, +1.324720240514420e+3,
                        +1.617548476343347e+4, +1.112729040439685e+2,
                        +1.074624783191125e+2, +8.835727765158191e+1,
                        +9.354078136054179e+1, +9.418142823531573e+1,
                        +1.040012390717851e+2, +6.861882624343235e+1,
                        +8.766654491283722e+1, +1.056007619389650e+2,
                        +7.738987569039419e+1, +1.041366366475571e+2])
    alpha_i = np.array([-6.743912502859256e+2, -3.973203432721332e+2,
                        -2.041233768918671e+3, -1.652917287299683e+3,
                        -1.783617639907328e+4, -5.887068595142284e+4,
                        -9.953255345514560e+3, -1.427131226068449e+3,
                        -3.256885197214938e+6, -2.924284515884309e+4,
                        -1.121774011188224e+3, -6.370088443140973e+4,
                        -1.008798413156542e+6, -8.837109731680418e+1,
                        -1.457246116408180e+2, -6.388286188419360e+1,
                        -2.195424319460237e+2, -6.719055740098035e+2,
                        -1.693747595553868e+2, -1.177598523430493e+1,
                        -4.596464999363902e+3, -1.738294585524067e+3,
                        -4.311715386228984e+1, -2.777743732451969e+2])
    alpha = np.array(alpha_r + alpha_i * 1j, dtype=np.complex128)
    n = A.shape[0]

    alpha0 = 2.258038182743983e-47

    k = 24

    y = np.array(n0, dtype=np.float64)
    for l in range(k):
        y = 2.0*np.real(alpha[l]*sla.spsolve(A*dt - theta[l]*sp.eye(n), y)) + y

    y *= alpha0

    Avogadro_Constant = 6.022140857E+23
    Electron_Coulomb = 1.6021766208E-19

    # 放射性活度，单位 Ci
    y_act = y * Avogadro_Constant * Coeff["Lambda"] / (3.7e10)

    
    # 衰变热，单位 Watts   
    y_q = y * Avogadro_Constant * Coeff["Lambda"] * Coeff["DecayHeat"] * Electron_Coulomb * 1.0e6

    # 放射性毒性，单位 Sv
    y_toxi = y * Avogadro_Constant * Coeff["Lambda"] * Coeff["Toxicity"]

    return y, y_act, y_q, y_toxi

