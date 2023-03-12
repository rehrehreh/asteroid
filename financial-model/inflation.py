# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 11:24:00 2023

@author: HelloWorld
"""

import pandas as pd
import cpi

startYear = 1994
endYear = 2023
startingAmount = 100

out = list()
for year in range(0,endYear-startYear):
    actYear = year+startYear
    actValue = cpi.inflate(startingAmount, startYear, to=actYear)
    out.append([actYear, actValue])
    
df=pd.DataFrame(out, columns=['year', 'value'])