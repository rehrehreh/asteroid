# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 17:02:25 2023

@author: mrehberg
"""
import math
import pandas as pd
import matplotlib.pyplot as plt
 
# ISP = 350
# v_e = ISP * 9.81 / 1000 
# m_dry=500
# dV1=2
# dV2=2
# m_asteroid=1000
# m_guess_1=100
# m_guess_2=100
# # This is the amount of propellant harvested from the asteroid
# # 1 corresponds to 100% propellant of the return prop comes from the asteroid
# Ap=1



def calculatePropellant(v_e,m_dry,dV1,dV2,m_asteroid,Ap,m_guess_1, m_guess_2):
    m_prop_1 = m_guess_1
    m_prop_2 = m_guess_2
    
    dV1_calc = v_e * math.log((m_dry + m_prop_1 + (1-Ap)*m_prop_2)/(m_dry + (1-Ap)*m_prop_2))
    dV2_calc = v_e * math.log((m_dry + m_prop_2 + m_asteroid)/(m_dry + m_asteroid))
    
    diff_dV1 = dV1 - dV1_calc
    diff_dV2 = dV2 - dV2_calc
    return m_prop_1, m_prop_2, diff_dV1, diff_dV2



# # 
# max_iter = 1000
# results = []
# for Ap in range(0,100):
#     Ap = Ap/100
#     for i in range(max_iter):
#         m_prop_1, m_prop_2, diff_dV1, diff_dV2 = calculatePropellant(v_e,m_dry,dV1,dV2,m_asteroid,Ap,m_guess_1, m_guess_2)
#         if abs(diff_dV1) +abs(diff_dV2) > 0.0001:
#             m_guess_1 = m_guess_1*(1+diff_dV1)
#             m_guess_2 = m_guess_2*(1+diff_dV2)
#         else: 
#             results.append((Ap, m_prop_1, m_prop_2))
#             break


# df = pd.DataFrame(results, columns=['Ap', 'propToAsteroid', 'propToEML'])
# df['launchMass'] = df['propToAsteroid'] + df['propToEML'] + m_dry 


# # Graphing # 

# fig, ax = plt.subplots(figsize=(12,6))
# ax.plot(df['Ap']*100, df['propToAsteroid'], label='propToAsteroid')
# ax.plot(df['Ap']*100, df['propToEML'], label='propToEML1')
# plt.axhline(y=m_asteroid, color='black', linestyle='--', label='massAsteroid')
# ax.set_xlabel('Percent of Return propellant harvested from asteroid (%)')
# ax.set_ylabel('Propellant Needed (kg)')
# ax.legend()
# ax.set_title(f'Dry Mass: {m_dry} kg, Asteroid Return: {m_asteroid} kg, dV2Ast: {dV1} km/s, dV2EML1: {dV2} km/s, ISP: {ISP} s')

