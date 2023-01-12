import pandas as pd
import numpy as np

global case
case = 'case1'

def readInputs(case):
    # read input file
    df_input=pd.read_csv(case+'/inputs.csv')
    # filter for names that start with a # - these are comments
    df_input=df_input[df_input['Name'].str[0]!='#']
    # fill NaN of numeric columns
    numericCols=['min', 'p10','p50','p90','max','numberDecimals']
    df_input[numericCols]=df_input[numericCols].fillna(-1)
    # Assign data types
    dataTypeDict={'numberDecimals': int,
                'min': float,
                'p10': float,
                'p50': float,
                'p90': float,
                'max': float}
    df_input= df_input.astype(dataTypeDict)

    return df_input


df_input=readInputs(case)
