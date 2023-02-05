# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 17:02:25 2023

@author: mrehberg
"""
import math
import pandas as pd
import matplotlib.pyplot as plt
 



def calculatePropellant(v_e,m_dry,dV1,dV2,m_asteroid,Ap,m_guess_1, m_guess_2):
    m_prop_1 = m_guess_1
    m_prop_2 = m_guess_2
    
    dV1_calc = v_e * math.log((m_dry + m_prop_1 + (1-Ap)*m_prop_2)/(m_dry + (1-Ap)*m_prop_2))
    dV2_calc = v_e * math.log((m_dry + m_prop_2 + m_asteroid)/(m_dry + m_asteroid))
    
    diff_dV1 = dV1 - dV1_calc
    diff_dV2 = dV2 - dV2_calc
    return m_prop_1, m_prop_2, diff_dV1, diff_dV2





# ISP = 420
# v_e = ISP * 9.81 / 1000 
# m_dry=900
# dV1=2.5
# dV2=2.5
# m_asteroid=1000
# m_guess_1=100
# m_guess_2=100
# # This is the amount of propellant harvested from the asteroid
# # 1 corresponds to 100% propellant of the return prop comes from the asteroid
# Ap=1


# # # 
# max_iter = 1000
# results = []
# for m in range(200,1200):
#     for i in range(max_iter):
#         m_prop_1, m_prop_2, diff_dV1, diff_dV2 = calculatePropellant(v_e,m,dV1,dV2,m_asteroid,Ap,m_guess_1, m_guess_2)
#         if abs(diff_dV1) +abs(diff_dV2) > 0.0001:
#             m_guess_1 = m_guess_1*(1+diff_dV1)
#             m_guess_2 = m_guess_2*(1+diff_dV2)
#         else: 
#             results.append((m, m_prop_1, m_prop_2))
#             break


# df = pd.DataFrame(results, columns=['dryMass', 'propToAsteroid', 'propToEML'])
# df['launchedPropellant'] = df['propToAsteroid'] + df['propToEML']*(1-Ap)
# df['NetPropellant'] = m_asteroid - df['launchedPropellant']

# # Graphing # 

# fig, ax = plt.subplots(figsize=(12,8))
# plt.rcParams.update({'font.size': 14})
# # ax.plot(df['Ap']*100, df['propToAsteroid'], label='propToAsteroid')
# # ax.plot(df['Ap']*100, df['propToEML'], label='propToEML1')
# ax.plot(df['dryMass'], df['NetPropellant'], label='NetPropellant', linewidth=5)
# plt.axhline(y=0, color='black', linestyle='--', label='AsteroidMaterialReturned')
# ax.set_xlabel('Dry Mass (kg)')
# ax.set_ylabel('Net Propellant [Harvested - Launched] (kg)')
# ax.legend()
# ax.set_title(f'Ap: {Ap*100}%, Asteroid Return: {m_asteroid} kg, dV2Ast: {dV1} km/s, dV2EML1: {dV2} km/s, ISP: {ISP} s')

