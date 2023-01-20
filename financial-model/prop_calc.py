# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 17:02:25 2023

@author: mrehberg
"""
import math 
v_e=2.942
m_dry=500
dV1=2
dV2=2.5
m_asteroid=1000
m_guess_1=100
m_guess_2=100
Ap=1
max_iter = 1000
def calculatePropellant(v_e,m_dry,dV1,dV2,m_asteroid,Ap,m_guess_1, m_guess_2):
    m_prop_1 = m_guess_1
    m_prop_2 = m_guess_2
    
    dV1_calc = v_e * math.log((m_dry + m_prop_1 + (1-Ap)*m_prop_2)/(m_dry + (1-Ap)*m_prop_2))
    dV2_calc = v_e * math.log((m_dry + m_prop_2 + m_asteroid)/(m_dry + m_asteroid))
    
    diff_dV1 = dV1 - dV1_calc
    diff_dV2 = dV2 - dV2_calc
    

    return m_prop_1, m_prop_2, diff_dV1, diff_dV2

for i in range(max_iter):
    m_prop_1, m_prop_2, diff_dV1, diff_dV2 = calculatePropellant(v_e,m_dry,dV1,dV2,m_asteroid,Ap,m_guess_1, m_guess_2)
    if abs(diff_dV1) +abs(diff_dV2) > 0.0001:
        m_guess_1 = m_guess_1*(1+diff_dV1)
        m_guess_2 = m_guess_2*(1+diff_dV2)
    else: 
        break

        