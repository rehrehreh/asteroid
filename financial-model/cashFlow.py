# Running cashFlow Model
import pandas as pd
import numpy as np
import random
from scipy.stats import norm
import matplotlib.pyplot as plt

global case, numYears
case = 'case1'
numYears = 25
numTrials = 1000

def readInputs(case):
    # read input file
    input=pd.read_csv(case+'/inputs.csv')
    # filter for names that start with a # - these are comments
    input=input[input['Name'].str[0]!='#'].reset_index(drop=True)
    # fill NaN of numeric columns
    numericCols=['varMin', 'p10','p50','p90','varMax','numberDecimals']
    input[numericCols]=input[numericCols].fillna(-1)
    # Assign data types
    dataTypeDict={'numberDecimals': int,
                'varMin': float,
                'p10': float,
                'p50': float,
                'p90': float,
                'varMax': float}
    input= input.astype(dataTypeDict)
    return input

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

def defineInputData(input):
    #check for duplicate variable names
    if max(input['varName'].value_counts())>1:
        raise ValueError("There are duplicate varNames")
    # Create an input dataframe for number of years specified
    var_dict=dict()
    for row in range(len(input)):
        use=input.loc[row,'use']
        varType = input.loc[row,'varType']
        varName=input.loc[row, 'varName']
        
        # Generate a random variable within the bounds
        if varType == 'numerical':
            if use == 'randomize':
                value = weibullDistribution(input.loc[row, ['p10', 'p50', 'p90', 'varMin', 'varMax']])
            else:
                value = input.loc[row, use]
            value = round(value, input.loc[row, 'numberDecimals'])
        
        # generate a random value from the list provided
        elif varType == 'categorical':
            if use == 'randomize':
                randomList = input.loc[row,'categorical'].split(',')
                value = random.choice(randomList)
            else:
                value = input.loc[row, use]
                
        var_dict.update({varName: value})
    return var_dict

def excavationSystemMass(var):
    # this function calculates the mass of the excavation system
    # the size is related to the number of scoops, water percent, and the water goal
    regolithNeeded = var['waterGoal']/var['waterPercent']
    excavationVolume = regolithNeeded/var['excavationNumScoops']/var['asteroidDensity']
    excavationSAtoVratio = var['excavationVolumeFactor']/excavationVolume**(1/3)
    excavationSystemSA = excavationSAtoVratio*excavationVolume
    excavationMass = var['excavationMaterialDensity'] * var['excavationMaterialThickness']*\
        excavationSystemSA * var['excavationMassFactor']   
    var.update({'excavationMass':round(excavationMass,0)})
    var.update({'regolithNeeded':round(regolithNeeded,0)})
    return 


inputs = readInputs(case)
outputs = pd.DataFrame()
for trial in range(numTrials):
    var = defineInputData(inputs)
    excavationSystemMass(var)
    outputs = pd.concat([outputs,pd.DataFrame.from_dict(var, orient='index').T])


### Graphing ###m
# redo everything with the mean variable 
# with the p10, mean, and p90 of the input variable
# plot these as bar graphs, sorted from largest to least deviated
# def tornado(output, outputVar, inputVars):
#     axs[0, 0].set_title('Output: {} vs {}'.format(outputVar,''.join(inputVars,', ')), fontsize=fs)
#     return

# tornado(outputs, 'excavationMass',  ['waterPercent', 'excavationMaterialDensity'])

outputs.plot.scatter('waterPercent','excavationMass')
    
    
    