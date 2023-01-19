# Running cashFlow Model
import pandas as pd
import numpy as np
import random
from scipy.stats import norm
import matplotlib.pyplot as plt

global case, numYears
case = 'case1'
numYears = 25
numTrials = 10

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

def defineInputData(inputs, force=None, forceVar = None):
    #check for duplicate variable names
    if max(inputs['varName'].value_counts())>1:
        raise ValueError("There are duplicate varNames")
    # Create an input dataframe for number of years specified
    var_dict=dict()
    for row in range(len(inputs)):
        varType = inputs.loc[row,'varType']
        varName=inputs.loc[row, 'varName']
        if force != None:
            if forceVar == varName:
                use = force
            else: 
                use = 'p50'
        else:
            use=inputs.loc[row,'use']
        # Generate a random variable within the bounds
        if varType == 'numerical':
            if use == 'randomize':
                value = weibullDistribution(inputs.loc[row, ['p10', 'p50', 'p90', 'varMin', 'varMax']])
            else:
                value = inputs.loc[row, use]
            value = round(value, inputs.loc[row, 'numberDecimals'])
        
        # generate a random value from the list provided
        elif varType == 'categorical':
            if use == 'randomize':
                randomList = inputs.loc[row,'categorical'].split(',')
                value = random.choice(randomList)
            else:
                value = inputs.loc[row, use]
                
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
    
    var.update({'excavationVolume':round(excavationVolume,2)})
    var.update({'excavationMass':round(excavationMass,0)})
    var.update({'regolithNeeded':round(regolithNeeded,0)})
    return 

def processingSystemMass(var):
    # This function is used to calculate the amount of energy needed to convert the asteroid material to water
    # It will then calculate the mass of the energy source and the processing equipment needed
    thermalEnergyPerKg = var['asteroidCp'] * (var['processEndTemp'] - var['processStartTemp']) / 1000 # kj/kg
    phaseChangePerKg = var['processEnthalpy'] / var['asteroidMW'] * 1000 #kJ/kg
    totalEnergyPerKg = thermalEnergyPerKg + phaseChangePerKg #kJ/kg
    totalEnergyPerBatch = totalEnergyPerKg * var['regolithNeeded'] / var['excavationNumScoops'] #kJ
    powerPerBatch = totalEnergyPerBatch / var['processTime'] * (1/24) * (1/60) * (1/60) #kW
    totalProcessingTime = var['processTime'] * var['excavationNumScoops'] # days

    # Mass of energy supply
    solarThermalRatio = 1 - var['processSolarPanelHeatRatio']
    totalSolarThermalMass = powerPerBatch * solarThermalRatio / var['solarThermalEnergyDensity'] / var['solarThermalEfficiency'] / (1 / var['asteroidMaxSunDistance']**2)
    totalSolarPanelMass = powerPerBatch * var['processSolarPanelHeatRatio'] / var['solarPanelEnergyDensity'] / var['solarPanelEfficiency'] / (1 / var['asteroidMaxSunDistance']**2)
    totalPowerMass = totalSolarPanelMass + totalSolarThermalMass

    # Mass of processing container
    # Assuming the batch size is the same for both systems
    processSAtoVRatio = var['processVolumeFactor']/var['excavationVolume']**(1/3)
    processSystemSA = processSAtoVRatio * var['excavationVolume']
    processContainerMass = var['processMaterialDensity'] * var['processMaterialThickness'] * processSystemSA * var['processMassFactor']
    
    totalProcessingMass = totalPowerMass + processContainerMass
    var.update({'totalProcessingMass':round(totalProcessingMass,2)})
    var.update({'powerPerBatch':round(powerPerBatch,3)})
    var.update({'totalPowerMass':round(totalPowerMass,2)})
    var.update({'totalEnergyPerBatch':round(totalEnergyPerBatch,2)})
    var.update({'totalProcessingTime':round(totalProcessingTime,2)})
    var.update({'processContainerMass':round(processContainerMass,2)})
    var.update({'totalEnergyPerKg':round(totalEnergyPerKg,2)})
    return

def summationVariables(var):
    totalPayloadMass = var['totalProcessingMass'] + var['excavationMass']
    var.update({'totalPayloadMass':round(totalPayloadMass,0)})
    return

def runSim(var):
    excavationSystemMass(var)
    processingSystemMass(var)
    summationVariables(var)
    return




### Graphing ###
def tornado(outputVar, inputs):
    #remove categorical variables for now
    df_inputs = inputs[inputs['varType']=='numerical'].reset_index(drop=True)
    # create a dataframe to house the sensitivity values
    senseCols = ['varMin', 'p10', 'p50', 'p90', 'varMax']
    sensitivities = pd.DataFrame(columns=['varName'] + senseCols)
    for row in range(len(df_inputs)):
        # run though all of the calculations for the min, p10, p50, p90, and max value
        varName = df_inputs.loc[row, 'varName']
        sensitivities.loc[row,'varName'] =  varName
        for sense in senseCols:
            #copy running statements from the trial loop
            # could put these in a for loop for single running eventually
            varTest = defineInputData(inputs, force=sense, forceVar=varName)
            runSim(varTest)
            sensitivities.loc[row,sense] = float(varTest[outputVar])
    
    sensitivities['sum'] = 0
    for col in senseCols:
        if col not in ['varName','p50']:
            sensitivities[col] = sensitivities[col] - sensitivities['p50']
            sensitivities['sum']+=abs(sensitivities[col])
    

    # ##ploting
    sensitivities = sensitivities.sort_values(by=['varMax'])
    sensitivities = sensitivities[sensitivities['sum']!=0].reset_index(drop=True)
    
    fig, ax = plt.subplots(figsize=(10,0.75*len(sensitivities)))
    ax.barh(sensitivities['varName'], sensitivities['varMin'] , align='center', color='firebrick', label='varMin')
    ax.barh(sensitivities['varName'], sensitivities['p10'] , align='center', color='palevioletred', label='p10')
    ax.barh(sensitivities['varName'], sensitivities['varMax'] , align='center', color='limegreen',label='varMax')
    ax.barh(sensitivities['varName'], sensitivities['p90'] , align='center', color='yellowgreen', label='p90')
    plt.axvline(x=0, color='black', linestyle='--')
    ax.set_yticks(sensitivities['varName'].to_list())
    plt.legend()
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_xlabel(f'Change to {outputVar}')
    meanVal = str(sensitivities['p50'].mean())
    ax.set_title(f'Mean {outputVar} is {meanVal}')
    return






inputs = readInputs(case)
outputs = pd.DataFrame()
for trial in range(numTrials):
    var = defineInputData(inputs)
    runSim(var)
    outputs = pd.concat([outputs,pd.DataFrame.from_dict(var, orient='index').T])
    
tornado('excavationMass', inputs)
tornado('totalEnergyPerKg', inputs)  
tornado('powerPerBatch', inputs)
tornado('totalEnergyPerKg', inputs)    
tornado('totalPayloadMass', inputs) 