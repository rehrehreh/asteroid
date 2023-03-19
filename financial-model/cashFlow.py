# Running cashFlow Model
import pandas as pd
import random
from scipy.stats import norm
import cpi
import matplotlib.pyplot as plt
from prop_calc import calculatePropellant
import os

global cases, saveFolder, cwd


def readInputs(cases):
    # read common inputs
    common = pd.read_csv('common-inputs.csv')
    common = processInputs(common)
    # read input file
    combine = common.copy()
    for case in cases:
        input=pd.read_csv(case+'/inputs.csv')
        input = processInputs(input)
        combine = pd.concat([combine, input])
    # combine these DFs, and keep the specific inputs
    combine = combine.drop_duplicates(subset=['varName'], keep='last').reset_index(drop=True)
    return combine

def processInputs(input):
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
    if 'asteroidWaterNeeded' in var:
        regolithNeeded = var['asteroidWaterNeeded']/var['waterPercent']
    else:
        regolithNeeded = var['waterGoal']/var['waterPercent']
    excavationVolume = regolithNeeded/var['excavationNumScoops']/var['asteroidDensity']
    excavationSAtoVratio = var['excavationVolumeFactor']/excavationVolume**(1/3)
    var['excavationSystemSA'] = excavationSAtoVratio*excavationVolume
    excavationMass = var['excavationMaterialDensity'] * var['excavationMaterialThickness'] * var['excavationSystemSA'] * var['excavationMassFactor']+ var['excvationMassA0']
    totalExcavationTime = var['excavationNumScoops'] * var['excavationTime']

    var.update({'excavationVolume':round(excavationVolume,2)})
    var.update({'excavationMass':round(excavationMass,0)})
    var.update({'regolithNeeded':round(regolithNeeded,0)})
    var.update({'totalExcavationTime':round(totalExcavationTime,2)})
    return 

def processingSystemMass(var):
    # This function is used to calculate the amount of energy needed to convert the asteroid material to water
    # It will then calculate the mass of the energy source and the processing equipment needed
    thermalEnergyPerKg = var['asteroidCp'] * (var['processEndTemp'] - var['processStartTemp']) / 1000 # kj/kg
    var['waterPerBatch'] = var['asteroidWaterNeeded'] / var['excavationNumScoops'] /var['processNumBatch']
    var['reactionEnergy'] = var['processEnthalpy'] * (var['waterPerBatch'] / .01802)  #kJ
    totalEnergyPerBatch = thermalEnergyPerKg * var['regolithNeeded'] / var['excavationNumScoops'] + var['reactionEnergy'] #kJ
    powerPerBatch = totalEnergyPerBatch / var['processTime'] * (1/24) * (1/60) * (1/60) #kW
    totalProcessingTime = var['processTime'] * var['excavationNumScoops'] * var['processNumBatch'] # days

    # Mass of energy supply
    solarThermalRatio = 1 - var['processSolarPanelHeatRatio']
    totalSolarThermalMass = powerPerBatch * solarThermalRatio / var['solarThermalEnergyDensity'] / var['solarThermalEfficiency'] / (1 / var['asteroidMaxSunDistance']**2)
    totalSolarPanelMass = powerPerBatch * var['processSolarPanelHeatRatio'] / var['solarPanelEnergyDensity'] / var['solarPanelEfficiency'] / (1 / var['asteroidMaxSunDistance']**2) / (1-var['solarAnnualDegradation']/100)**(var['designYears'])
    totalPowerMass = totalSolarPanelMass + totalSolarThermalMass

    # Mass of processing container
    # Assuming the batch size of the processing system is on a daily basis
    var['processingVolume'] = var['excavationVolume'] / var['processNumBatch'] / var['processTime']
    processSAtoVRatio = var['processVolumeFactor']/var['processingVolume']**(1/3)
    var['processSystemSA'] = processSAtoVRatio * var['processingVolume']
    processContainerMass = var['processMaterialDensity'] * var['processMaterialThickness'] * var['processSystemSA'] * var['processMassFactor']
    
    totalProcessingMass = totalPowerMass + processContainerMass
    var.update({'totalProcessingMass':round(totalProcessingMass,2)})
    var.update({'powerPerBatch':round(powerPerBatch,3)})
    var.update({'totalPowerMass':round(totalPowerMass,2)})
    var.update({'totalEnergyPerBatch':round(totalEnergyPerBatch,2)})
    var.update({'totalProcessingTime':round(totalProcessingTime,2)})
    var.update({'processContainerMass':round(processContainerMass,2)})
    return

def spacecraftPowerCalculation(var):
    var['baseSolarPanelMass'] = (var['spacecraftBasePower'] + var['cryocoolerPower'] * var['cryocooler']) \
                                 * var['solarPanelEnergyDensity'] / (1 / var['asteroidMaxSunDistance']**2) / (1-var['solarAnnualDegradation']/100)**(var['designYears'])
    var['batteryMass'] = var['spacecraftBasePower'] * 1000 *var['batteryTime'] / var['batteryDensity']
    
    ## Cost RD
    var['costPower1Rd'] = var['cerEPSA1'] * (var['batteryMass'] + var['baseSolarPanelMass'])
    var['costPower2Rd'] = var['cerEPSA2'] * ((var['batteryMass'] + var['baseSolarPanelMass']) * var['spacecraftBasePower'])**var['cerEPSE1']
    var['costPowerRd'] = var['costPower1Rd'] + var['costPower2Rd']
    
    ## Cost TFU
    var['costPowerTfu'] = var['cerTfuEPSA1'] * (var['batteryMass'] + var['baseSolarPanelMass'])**var['cerTfuEPSE1']
    return

def spacecraftThermalCalculation(var):
    var['spacecraftSurfaceArea'] = (var['processSystemSA'] + var['excavationSystemSA'] + var['tankSA']) * var['spacecraftSAFactor']
    var['massMLI'] = var['spacecraftSurfaceArea'] * var['mliDensity']


    var['thermalLoad'] = var['radiationPower']
    var['radiatorArea'] = var['thermalLoad'] * 1000 / (5.67e-8 * var['radiationRejectionTemperature']**4)
    var['massRadiator'] = var['radiatorArea'] * var['radiatorDensity']
    
    ## Cost
    
    return

def spacecraftCommandAndDataHandlingCalculation(var):
    var['massCDH'] = var['cdhMass']
    
    ## Cost

    return

def spacecraftStructureCalculation(var):
    ## Mass is a percent of dry mass
    # this mass is not added to the dry mass, it is just for cost calculations
    var['structureMass'] = var['dryMass']*var['structureDryMassFactor']
    ## Cost
    var['costStructureRd'] = var['cerStructureA1'] * var['structureMass']**var['cerStructureE1']
    return

def spacecraftGuidanceAndNavigationCalculation(var):
    var['massGNC'] = var['gncMass']
    
    ## Cost
    
    return

def spacecraftCommunicationsCalculation(var):
    var['massComms'] = var['commMass']
    return

def spacecraftAttitudeControlCalculation(var):
    var['massAtt'] = var['attMass']
    return

def spacecraftElectroysis(var):
    if var['electrolysis'] > 0:
        var['hydrogenProductionMass'] = var['massPropToEML1'] / 18.02 * 2.02
        var['elecEnergy'] = var['hydrogenProductionMass'] * var['elecHydrogenEnergy']
        if var['elecPower'] < (var['powerPerBatch'] * var['processSolarPanelHeatRatio']):
            var['powerElec'] = var['powerPerBatch'] * var['processSolarPanelHeatRatio']
            var['elecPowerMass'] = 0
        else:
            var['powerElec'] = var['elecPower']
            var['elecPowerMass'] = var['powerElec']  / var['solarPanelEnergyDensity'] / var['solarPanelEfficiency'] / (1 / var['asteroidMaxSunDistance']**2) / (1-var['solarAnnualDegradation']/100)**(var['designYears'])
        var['timeElectrolysis'] = var['elecEnergy'] / var['powerElec']  /24
        var['massElectrolysis'] = var['elecPowerMass'] + var['elecMass']
    else:
        var['timeElectrolysis'] = 0
        var['massElectrolysis'] = 0
        var['powerElec'] = 0
        var['elecPowerMass'] = 0
    return

def spacecraftPropulsionCalulation(var):
    if 'totalPropellant' in var:
        totalPropellant = var['totalPropellant']
    else:
        totalPropellant = 200 #kg
    var['propVolume'] = totalPropellant / var['propDensity'] / var['tankNumber']
    var['tankSA'] = var['tankVolumeFactor'] / var['propVolume']**(1/3) * var['tankNumber']
    var['massTank'] = var['tankPressurantMassFactor'] * (var['tankSA'] * var['tankMaterialThickness'] * var['propTankDensity'])
    var['massEngine'] = var['engineMass']
    return

def cryocooler(var):
    if var['cryocooler'] == 0:
        var['cryocoolerMass'] = 0
        var['cryocoolerPower'] = 0
        var['cryocoolerRDCost'] = 0
        var['cryocoolerTFUCost'] = 0
    return 

def waterFiltration(var):
    if var['waterFilter'] == 0:
        var['waterFilterMass'] = 0
        var['waterFilterPower'] = 0
        var['waterFilterLoss'] = 0
    
    if var['waterFilterPower'] < (var['powerPerBatch'] * var['processSolarPanelHeatRatio']):
        var['waterFilterPowerMass'] = 0
    else:
        var['waterFilterPowerMass'] = var['waterFilterPower']  / var['solarPanelEnergyDensity'] / var['solarPanelEfficiency'] / (1 / var['asteroidMaxSunDistance']**2) / (1-var['solarAnnualDegradation']/100)**(var['designYears'])
    
    var['massWaterFilter'] = var['waterFilterPowerMass'] + var['waterFilterMass']
    # water filter loss accounted for in calculateTotalPropellant
    return
 
def summationVariables(var):
    var['totalPayloadMass'] = var['totalProcessingMass'] + var['excavationMass']
    var['totalStayDays'] = var['totalProcessingTime'] + var['excavationTime']
    var['dryMass'] = var['totalProcessingMass'] + var['excavationMass'] + var['baseSolarPanelMass'] +\
        var['cdhMass'] +  var['massRadiator'] + var['massMLI'] + var['massGNC'] + var['massEngine'] +\
            var['massTank'] + var['commMass'] + var['massAtt'] + var['batteryMass'] + var['massElectrolysis'] +\
                var['cryocoolerMass'] + var['massWaterFilter']
    var['dryMass']  = var['dryMass'] * (1 + var['structureDryMassFactor'] + var['dryMassMargin'])
    var['netProp'] = var['waterGoal'] - var['massPropToAsteroid']
    var['launchMass'] = var['dryMass'] + var['launchedPropellant']
    var['percentPropSold'] = var['netProp'] / var['waterGoal']
    return

def simOrder(var):
    calculateTotalPropellant(var)
    excavationSystemMass(var)
    processingSystemMass(var)
    spacecraftPowerCalculation(var)
    spacecraftPropulsionCalulation(var)
    spacecraftGuidanceAndNavigationCalculation(var)
    spacecraftThermalCalculation(var)
    spacecraftAttitudeControlCalculation(var)
    spacecraftCommunicationsCalculation(var)
    spacecraftCommandAndDataHandlingCalculation(var)
    spacecraftElectroysis(var)
    cryocooler(var)
    waterFiltration(var)
    summationVariables(var)
    spacecraftStructureCalculation(var)
    return




def calculateTotalPropellant(var):
    isp = var['engineISP']
    dV1 = var['deltaVtoAsteroid']
    dV2 = var['deltaVtoEML1']
    v_e = isp * 9.81 / 1000 
    if 'dryMass' not in var:
         dryMass = var['dryMassGuess'] * var['spacecraftDryMassFactor'] * (1 + var['dryMassMargin'])
    else: 
        dryMass = var['dryMass']
    
    Ap = var['asteroidPropellantRatio']
    massProp1Guess = 1
    massProp2Guess = 2
    massAsteroid = var['waterGoal']
    maxIterations = 100
    convergence = 0.01
    for i in range(maxIterations):
        massProp1, massProp2, diff_dV1, diff_dV2 = calculatePropellant(v_e, dryMass, dV1, dV2, massAsteroid, Ap, massProp1Guess, massProp2Guess)
        if abs(diff_dV1) +abs(diff_dV2) > convergence:
            massProp1Guess = massProp1Guess*(1+diff_dV1)
            massProp2Guess = massProp2Guess*(1+diff_dV2)
        else: 
            break
    var['totalPropellant'] = massProp1 + massProp2 
    var['launchedPropellant'] = massProp1 + (1- var['asteroidPropellantRatio']) * massProp2
    var['asteroidWaterNeeded'] = var['waterGoal'] + massProp2 *(var['asteroidPropellantRatio']) + var['waterFilterLoss']*var['waterFilter']
    
    var['netPropellant'] = var['waterGoal'] - var['launchedPropellant']
    
    var.update({'massPropToAsteroid':round(massProp1,0)})
    var.update({'massPropToEML1':round(massProp2,0)})
    return


def calculateSpacecraftCost(var):

    var['rdte_costCer_spacecraft'] = var['dryMass'] * var['cerDryMassA1'] + var['costStructureRd'] + var['costPowerRd'] + var['cryocoolerRDCost']
    var['rdte_costCer_iatX1'] = var['rdte_costCer_spacecraft'] 
    var['rdte_costCer_iat'] = var['cerIATA0'] + var['cerIATA1'] * var['rdte_costCer_iatX1']
    # same X1 as IAT, integration assembly and test
    var['rdte_costCer_programLevelCostCER'] = var['cerProgramA1'] * var['rdte_costCer_iatX1'] **(var['cerProgramE1'])
    
    var['rdte_costCer_totalCost'] = var['rdte_costCer_programLevelCostCER'] + var['rdte_costCer_iat'] + var['rdte_costCer_spacecraft']
    
    #Theoretical First Unit Cost
    var['tfu_costCer_spacecraft'] = var['dryMass'] * var['cerTFUDryMassA1']
    var['tfu_costCer_iatX1'] = var['dryMass'] 
    var['tfu_costCer_iat'] = var['cerTFUIATA1'] * var['tfu_costCer_iatX1']
    
    var['tfu_costCer_programX1'] = var['tfu_costCer_spacecraft'] + var['tfu_costCer_iat']
    var['tfu_costCer_programLevelCostCER'] = var['cerTFUProgramA1'] * var['tfu_costCer_programX1']
    
    var['tfu_costCer_totalCost'] = var['tfu_costCer_programLevelCostCER'] + var['tfu_costCer_iat'] + var['tfu_costCer_spacecraft'] + var['costPowerTfu'] +var['cryocoolerTFUCost']
    var['tfu_costCer_totalCost'] = var['tfu_costCer_totalCost'] * var['cerTFUTRLFactor']
    
    var['totalCost'] = var['tfu_costCer_totalCost'] + var['rdte_costCer_totalCost']
    
    #convert to 2023 $
    var['rdte_costCer_totalCost'] = round(cpi.inflate(var['rdte_costCer_totalCost'],2000,to=2022),0)
    var['tfu_costCer_totalCost'] =round(cpi.inflate(var['tfu_costCer_totalCost'], 2000, to=2022),0)
    var['totalCost'] =round(cpi.inflate(var['totalCost'], 2000, to=2022) ,0)
    var['costPerKgWater'] = var['totalCost'] / var['waterGoal']
    return


def costMetrics(var):
    var['launchCost'] = (var['dryMass'] + var['launchedPropellant'])*var['launchPrice']
    var['costLaunchAndTFU'] = var['launchCost'] + var['tfu_costCer_totalCost']
    var['costLaunchTotal'] = var['totalCost'] + var['launchCost']
    var['soldProp'] = (var['roundTrips']-1)*var['netProp'] + var['waterGoal']
    var['recoverySalePrice'] = var['costLaunchAndTFU'] / var['soldProp']
    var['recoverySalePriceWithRD'] = var['costLaunchTotal'] / var['soldProp']
    return


    
def runSim(var):
    # these variables determine the maximum number of iterations and convergence minimum for the dry mass calculations
    maxIterations = 1000
    convergence = 1 #kg
    propGuess = 200
    for iter in range(maxIterations):
        if iter == 0:
            simOrder(var)
            propCheck = abs(propGuess - var['totalPropellant'])
        else:
            if propCheck > convergence:
                propGuess = var['totalPropellant']
                simOrder(var)
            else:
                break
    calculateSpacecraftCost(var)
    costMetrics(var)
    # 
    return


#####################################################################################################
#######################################################################################################
### Graphing ###
def addlabels(sensitivities, df_inputs, var):
    for i in range(len(sensitivities)):
        if i > 12:
            break
        y = sensitivities[var][i]
        varName = sensitivities['varName'][i]
        value = df_inputs[df_inputs['varName']==varName][var].values[0]
        plt.text(y, i, value, ha = 'center', fontsize=11, 
                 Bbox = dict(facecolor = 'white', alpha =.6))

def tornado(outputVar, inputs, maxEffect=0, saveFile=None, maxRound=None):
    plt.rcParams.update({'font.size': 16})
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
    sensitivities = sensitivities.sort_values(by=['sum'], ascending=False)
    sensitivities = sensitivities[sensitivities['sum']>=maxEffect].reset_index(drop=True)
    
    fig, ax = plt.subplots(figsize=(10,0.75*len(sensitivities)))
    ax.barh(sensitivities['varName'], sensitivities['varMin'] , align='center', color='firebrick', label='varMin')
    addlabels(sensitivities, df_inputs, 'varMin')
    ax.barh(sensitivities['varName'], sensitivities['p10'] , align='center', color='palevioletred', label='p10')
    addlabels(sensitivities, df_inputs, 'p10')
    ax.barh(sensitivities['varName'], sensitivities['varMax'] , align='center', color='limegreen',label='varMax')
    addlabels(sensitivities, df_inputs, 'varMax')
    ax.barh(sensitivities['varName'], sensitivities['p90'] , align='center', color='yellowgreen', label='p90')
    addlabels(sensitivities, df_inputs, 'p90')
    plt.axvline(x=0, color='black', linestyle='--')
    ax.set_yticks(sensitivities['varName'].to_list())
    plt.legend()
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_xlabel(f'Change to {outputVar}')
    medianVal = str(sensitivities['p50'].median())
    ax.set_title(f'Base {outputVar} is {medianVal}')
    if saveFile == True:
        plt.savefig(saveFolder + f'{outputVar}.png')
    return

def plottingOutputCorellations(outputs, x, y, ylim=None, xlim=None, saveFile=None):
    fig, ax = plt.subplots(figsize=(10,6))
    ax.scatter(outputs[x], outputs[y])
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    if ylim!=None:
        ax.set_ylim(ylim)
    if xlim!=None:
        ax.set_xlim(xlim)
    casesTitle = ', '.join(cases)
    ax.set_title(casesTitle)
    # ax.plot(np.unique(outputs[x]), np.poly1d(np.polyfit(outputs[x], outputs[y], 1))(np.unique(outputs[x])), color='r')
    if saveFile==True:
        plt.savefig(saveFolder + f'{x}_vs_{y}.png')
    return

def singleVariableRange(inputs,inputVar,varRangeMin,varRangeMax,factor=1):
    outputs = pd.DataFrame()
    for a in range(int(varRangeMin * factor), int(varRangeMax * factor)):
        a_in = a / factor
        var = defineInputData(inputs, force = 'Yes')
        var[inputVar] = a_in
        runSim(var)
        outputs = pd.concat([outputs,pd.DataFrame.from_dict(var, orient='index').T])    
    return outputs

def getpValues(outputs, outVars):
    pVals = pd.DataFrame()
    for outVar in outVars:
        temp = {}
        temp['var'] = outVar
        temp['min'] = min(outputs[outVar]).round(1)
        temp['p10'] = outputs[outVar].quantile(0.1).round(1)
        temp['p50'] = outputs[outVar].quantile(0.5).round(1)
        temp['p90'] = outputs[outVar].quantile(0.9).round(1)
        temp['max'] = max(outputs[outVar]).round(1)
        df_temp = pd.DataFrame.from_dict(temp, orient='index').T
        pVals = pd.concat([pVals, df_temp])
    return pVals

###################################################################################################################
###################################################################################################################



def main(inputs, running):
    ## Run a monte Carlo
    if running == 'MonteCarlo':
        outputs = pd.DataFrame()
        numTrials = 1000
        for trial in range(numTrials):
            var = defineInputData(inputs)
            runSim(var)
            #sanity check for results
            if var['dryMass']<6000:
                outputs = pd.concat([outputs,pd.DataFrame.from_dict(var, orient='index').T])
        outputs.to_csv(cwd + '\\outputs\\simData.csv', float_format='%.3f')
        
        # get the distributed outputs for Joes model
        outputs_to_Joe = ['dryMass', 'launchMass', 'netProp', 'totalStayDays', 'rdte_costCer_totalCost', 'tfu_costCer_totalCost']
        pVals = getpValues(outputs, outputs_to_Joe)
        print(pVals)
        pVals.to_csv(cwd +'\\outputs\\outputs_to_joe.csv', index=False)
        plottingOutputCorellations(outputs,  y='netProp', x='dryMass',  saveFile=True)
    
    ### run a asingle variable run
    if running == 'Range':
        minRange = 1000
        maxRange = 3000
        outputs = singleVariableRange(inputs,'waterGoal', minRange, maxRange, .1)
        plottingOutputCorellations(outputs, y='dryMass', x='waterGoal', xlim=[minRange,maxRange], saveFile=True)
        plottingOutputCorellations(outputs, y='netProp', x='waterGoal', xlim=[minRange,maxRange], saveFile=True)
        plottingOutputCorellations(outputs, y='percentPropSold', x='waterGoal', xlim=[minRange,maxRange], saveFile=True)
        plottingOutputCorellations(outputs, y='rdte_costCer_totalCost', x='waterGoal', xlim=[minRange,maxRange],  saveFile=True)
        plottingOutputCorellations(outputs, y='tfu_costCer_totalCost', x='waterGoal', xlim=[minRange,maxRange],  saveFile=True)
        plottingOutputCorellations(outputs, y='recoverySalePrice', x='waterGoal', xlim=[minRange,maxRange],  saveFile=True)
    
    
    
    ## Tornado Plots
    if running == 'Tornado':
        # tornado('dryMass', inputs, 20, saveFile=True)  
        #tornado('excavationMass', inputs, 5)
        #tornado('totalProcessingMass', inputs, 20)    
        # tornado('netProp', inputs, 50, saveFile=True)
        tornado('recoverySalePrice',inputs,10,saveFile=True)
        tornado('recoverySalePriceWithRD',inputs,40,saveFile=True)
        # tornado('tfu_costCer_totalCost', inputs, 5000, saveFile=True)
    return

if __name__ =='__main__':
    cwd = os.path.dirname(os.path.realpath(__file__))
    # Can specify as many as desired
    # But should only ever have a single engine and a single solar method
    # If you do not want any cases, specify 'base'
    cases = ['case-solar-thermal' ,'case-engine-dirtySteam']
    inputs = readInputs(cases)
    outputs = pd.DataFrame()
    
    
    # this is the type of analysis to run
    running = 'Tornado'
    saveFolder= cwd + f'\\figures\\{running}\\'
    main(inputs, running)