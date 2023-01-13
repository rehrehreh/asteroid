import pandas as pd
import numpy as np
import random
from scipy.stats import norm

global case
case = 'case1'

def readInputs(case):
    # read input file
    df_input=pd.read_csv(case+'/inputs.csv')
    # filter for names that start with a # - these are comments
    df_input=df_input[df_input['Name'].str[0]!='#']
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

def weibullDistribution(p10, p50, p90, varMin, varMax, number_of_samples):
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
        samples.append(sample)
    return samples

def defineInputData(df_input):
    

# varMin = 1100
# p10_in = 1350
# p50_in = 1450
# p90_in = 1500
# varMax = 2000
# numberofsamples = 10000
# data = weibullDistribution(p10_in, p50_in, p90_in, varMin, varMax, 1)
# pd.DataFrame(data).hist(bins=50)

# p10_out = np.percentile(data,10)
# p50_out = np.percentile(data,50)
# p90_out = np.percentile(data,90)




df_input=readInputs(case)
