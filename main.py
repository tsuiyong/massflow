#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
Author: Yong Cui @ DRP.SINAP
Date: 13:36 2018-9-12

Description
===========

This script is used to calculate the mass flow for fuel cycle code TMSR-NUDEAS.
"""

import datetime as dt
from func import alias2id
from func import mole2kg
from func import fission_product
from func import transuranium
from func import nuc1
from func import nuc2
from func import nuclide_evolution
import numpy as np
from func import paint_core
from func import paint_activity
from func import radio_activity
import pandas as pd
from func import required_nuclides
from func import paint_toxicity

start = dt.datetime.now()
df = pd.read_table('input.txt', sep='\s+')
power_density = df.iloc[0, 1]
volume = float(df.iloc[1, 1])
eta = float(df.iloc[2, 1])
Uranium_alias = df.iloc[3, 1]
year = float(df.iloc[4, 1])

enrichment = 1  # 该值后期应从数据库读入
power_electricity = float(power_density) * volume * eta / 1000

f = open("mass_flow.txt", "w")

if Uranium_alias == "U233" or Uranium_alias == "Pu239":
    f.write("uranium_enrichment(%):" + "\n")
    f.write(str(0) + "\n")
    f.write("enrichment_uranium_need(kg/GWe):" + "\n")
    f.write(str(0) + "\n")
else:
    f.write("uranium_enrichment(%):" + "\n")
    f.write(str(enrichment) + "\n")
    id = alias2id("U235")
    U235_mole = nuc1(power_density, id)
    U235_kg = mole2kg(U235_mole, id)
    id = alias2id("U238")
    U238_mole = nuc1(power_density, id)
    U238_kg = mole2kg(U238_mole, id)
    Uranium_need = U235_kg + U238_kg
    f.write("enrichment_uranium_need(kg/GWe):" + "\n")
    f.write(str(Uranium_need / power_electricity) + "\n")

nuclide_evol = nuclide_evolution(power_density, alias2id(Uranium_alias))
paint_core(year, Uranium_alias, power_density, mole2kg(nuclide_evol, alias2id(Uranium_alias)))

id = alias2id("Th232")
Thorium_mole = nuc1(power_density, id)
Thorium_need = mole2kg(Thorium_mole, id)
f.write("thorium_need(kg/GWe):" + "\n")
f.write(str(Thorium_need / power_electricity) + "\n")
nuclide_evol = nuclide_evolution(power_density, alias2id("Th232"))
paint_core(year, "Th232", power_density, mole2kg(nuclide_evol, alias2id("Th232")))

print("Calculating the mass of Fission Products(Kg)......")

mass_fission_products = fission_product(power_density)
f.write("mass_fission_products(kg/GWe):" + "\n")
f.write(str(mass_fission_products / power_electricity) + "\n")

isotope_kg = [mole2kg(nuc2(power_density, id), id) for id in range(1610, 1623)]
mass_depleted_uranium = np.sum(isotope_kg)
f.write("mass_depleted_uranium(kg/GWe):" + "\n")
f.write(str(mass_depleted_uranium / power_electricity) + "\n")

isotope_kg = [ mole2kg(nuc2(power_density, id), id) for id in range(1592, 1601)]
mass_depleted_thorium = np.sum(isotope_kg)
f.write("mass_depleted_thorium(kg/GWe):" + "\n")
f.write(str(mass_depleted_thorium / power_electricity) + "\n")

print("Calculating the mass of TRansUranium(Kg)......")

mass_TRU = transuranium(power_density)
f.write("mass_TRU(kg/GWe):" + "\n")
f.write(str(mass_TRU / power_electricity) + "\n")

print("Calculating the radio activity and TRU toxicity after the reactor is shut down......")
mass_off = [nuc2(power_density, i) for i in range(1693)]

decay_time = 100
total_act = radio_activity(decay_time, mass_off)
f.write("radioactive_100y_spent_fuel(Ci/GWe-yr):" + "\n")
f.write(str(total_act / power_electricity / year) + "\n")

paint_activity(decay_time, mass_off)
paint_toxicity(decay_time, mass_off)

decay_time = 1000000
total_act = radio_activity(decay_time, mass_off)
f.write("radioactive_10000thy_spent_fuel(Ci/GWe-yr):" + "\n")
f.write(str(total_act / power_electricity / year) + "\n")

paint_activity(decay_time, mass_off)
paint_toxicity(decay_time, mass_off)

f.close()

print("Calculating the mass(kg) and decay heat(W) of required nuclides......")

required_nuclides(mass_off)

end = dt.datetime.now()

print("Total computing time: " + str(end - start))
