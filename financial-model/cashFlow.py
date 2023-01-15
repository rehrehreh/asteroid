import pandas as pd
import numpy as np
import random
from scipy.stats import norm

global case, numYears
case = 'case1'
numYears = 25

def readInputs(case):
    # read input file
    df_input=pd.read_csv(case+'/inputs.csv')
    # filter for names that start with a # - these are comments
    df_input=df_input[df_input['Name'].str[0]!='#'].reset_index(drop=True)
    # fill NaN of numeric columns
    numericCols=['varMin', 'p10','p50','p90','varMax','numberDecimals']
    df_input[numericCols]=df_input[numericCols].fillna(-1)
    # Assign data types
    dataTypeDict={'numberDecimals': int,
                'varMin': float,
                'p10': float,
                'p50': float,
                'p90': float,
                'varMax': float}
    df_input= df_input.astype(dataTypeDict)
    return df_input

def weibullDistribution(input_list):
    [p10, p50, p90, varMin, varMax] = input_list
    number_of_samples = 1
    b_mark = ((float(p90) - float(p50)) / (float(p50) - float(p10)))
    samples = []
    for i in range(number_of_samples):
        rand_numb = random.random()
        factor = norm.ppf(rand_numb, 0, 0.780304146072379)
        if 0.9999 < b_mark < 1.0001: 
            sample = p50 + (p90 - p50) * factor
        else:
            sample = p50 + (p90 - p50)*((b_mark**factor - 1)/(b_mark - 1))
        sample=min(max(sample, varMin), varMax)
        #samples.append(sample)
    return sample

def defineInputData(df_input):
    #check for duplicate vavriable names
    if max(df_input['varName'].value_counts())>1:
        raise ValueError("There are duplicate varNames")
    # Create an input dataframe for number of years specified
    df=pd.DataFrame(range(numYears), columns=['year'])
    for row in range(len(df_input)):
        use=df_input.loc[row,'use']
        varName=df_input.loc[row, 'varName']
        # Generate a random variable within the bounds
        if use == 'randomize':
            value = weibullDistribution(df_input.loc[row, ['p10', 'p50', 'p90', 'varMin', 'varMax']])
        else:
            value = df_input.loc[row, use]
        value = round(value, df_input.loc[row, 'numberDecimals'])
        df[varName] = value
    return df



df_input=readInputs(case)
df = defineInputData(df_input)