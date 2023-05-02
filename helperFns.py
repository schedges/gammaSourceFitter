#Helper functions for the gamma source fitter
#We avoid any Roofit here out of concern for our observable not transferring over. 
#Helper functions include
#   Parser
#   Smearer
#   User guess interface
import re
import numpy as np
import ROOT
import array
import sys
import os

#Spits an error and the corresponding line if we can't parse it
def throwParserError(message,line):
    print(message)
    print(line)
    sys.exit(-1)

#Remove any formatting characters user inputted, and break into parameter name and arguments. Args stored in a list if more than one
def parseLine(line):
    #Remove formatting characters
    line=re.sub("\[| |\"|\]","",line)
    #Split based on colon
    lineParts=line.split(":")
    if not len(lineParts)==2:
        throwParserError("Tried to split line based on ':' but : sign not found",line)
    #Load into parameter name and arguments
    argName=lineParts[0]
    args=lineParts[1].split(",")
    if len(args)==0:
        throwParserError("No arguments found to right of : sign in line",line)
    return argName,args

#Checks if the arguments are valid (not missing, has correct # of values, and valid option
def checkIfInputValid(argName,arg,nValues,validOptions):
    #Catch missing input
    if arg=="":
        print("Error! No "+argName+" parameter found")
        sys.exit(-1)
    #Catch input w/fewer than required arguments
    if not isinstance(arg, list):
        length=1
    else:
        length=len(arg)
    if length<nValues:
        print("Error! "+argName+" requires at least "+str(nValues)+" arguments")
        print("You put: ")
        print(arg)
        sys.exit(-1)
    #If a list of validOptions is supplied, make sure the user argument is in that list
    if not validOptions==[]:
        if not arg in validOptions:
            print("Error! Valid "+argName+" options are:")
            print(validOptions)
            print("You put: ")
            print(arg)
            sys.exit(-1)

#Checks if the selected range is valid
def checkIfRangeValid(argName,arg):
    if len(arg)>1:
        if arg[0]>arg[1]:
            print("Error on parameter "+argName+", arguments not in ascending order")
            print(arg)
            sys.exit(-1) 

def checkIfFileExists(fname):
    if not os.path.exists(fname):
        print("Error: "+fname+" could not be found")
        sys.exit(-1)
        
def checkInputForErrors(argDict):
    #Make sure we have required arguments
    requiredArgs=["title","runType","binningType","binningParameter","fitRange","alphaRange",
                  "betaRange","gammaRange","slopeRange","offsetRange","simScalingRange","bgndName","bgndFilePath","bgndTreeName",
                  "bgndEnergyBranchName","bgndNormalization","bgndNormalizationParameter","sourceNames",
                  "sourcePaths","sourceTreeName","sourceEnergyBranchName","sourceNormalization",
                  "sourceNormalizationParameter","simNames","simPaths","simTreeName","simBranchName",
                  "nWalkers","nBurnInSteps","nSteps"]
    for arg in requiredArgs:
        if not arg in argDict:
            print("Missing required argument "+arg)
            print(argDict)
            print(sys.exit())

    #Check meta data
    #0. Title
    argDict["title"]=argDict["title"][0]
    checkIfInputValid("title",argDict["title"],1,[])
    #1. runType
    argDict["runType"]=argDict["runType"][0]
    checkIfInputValid("runType",argDict["runType"],1,["fit","guess"])
    #2. binningType   
    argDict["binningType"]=argDict["binningType"][0]
    checkIfInputValid("binningType",argDict["binningType"],1,["fixed","variable"])
    #3. binningParameter
    argDict["binningParameter"]=[float(i) for i in argDict["binningParameter"]]
    if argDict["binningType"]=="fixed":
        argDict["binningParameter"]=argDict["binningParameter"][0]
    else:
        print("")
        #TODO For variable binning, check range appropriate
    #5. fitRange
    argDict["fitRange"]=[float(i) for i in argDict["fitRange"]]
    checkIfInputValid("fitRange",argDict["fitRange"],2,[])
    checkIfRangeValid("fitRange",argDict["fitRange"])
    #6. alphaRange
    argDict["alphaRange"]=[float(i) for i in argDict["alphaRange"]]
    checkIfInputValid("alphaRange",argDict["alphaRange"],1,[])
    checkIfRangeValid("alphaRange",argDict["alphaRange"])    
    #7. betaRange
    argDict["betaRange"]=[float(i) for i in argDict["betaRange"]]
    checkIfInputValid("betaRange",argDict["betaRange"],1,[])
    checkIfRangeValid("betaRange",argDict["betaRange"])
    #8. alphaRgammaRangeange
    argDict["gammaRange"]=[float(i) for i in argDict["gammaRange"]]
    checkIfInputValid("gammaRange",argDict["gammaRange"],1,[])
    checkIfRangeValid("gammaRange",argDict["gammaRange"])
    #9. slopeRange
    argDict["slopeRange"]=[float(i) for i in argDict["slopeRange"]]
    checkIfInputValid("slopeRange",argDict["slopeRange"],1,[])
    checkIfRangeValid("gammaRaslopeRangenge",argDict["slopeRange"])
    #10. offsetRange
    argDict["offsetRange"]=[float(i) for i in argDict["offsetRange"]]
    checkIfInputValid("offsetRange",argDict["offsetRange"],1,[])
    checkIfRangeValid("offsetRange",argDict["offsetRange"])
    #11. simScalingRange
    argDict["simScalingRange"]=[float(i) for i in argDict["simScalingRange"]]
    checkIfInputValid("simScalingRange",argDict["simScalingRange"],1,[])
    checkIfRangeValid("simScalingRange",argDict["simScalingRange"])
    #12. maxSimEntries - optional
    if not "maxSimEntries" in argDict:   
        argDict["maxSimEntries"]=[-1]
    argDict["maxSimEntries"]=int(argDict["maxSimEntries"][0])
    #13. nWalkers
    argDict["nWalkers"]=int(argDict["nWalkers"][0])
    #13. nBurnInSteps
    argDict["nBurnInSteps"]=int(argDict["nBurnInSteps"][0])
    #13. nSteps
    argDict["nSteps"]=int(argDict["nSteps"][0])
    
    #Check Bgnd Data
    #1. bgndName
    argDict["bgndName"]=argDict["bgndName"][0]    
    #2. Bgnd file path
    argDict["bgndFilePath"]=argDict["bgndFilePath"][0]    
    checkIfFileExists(argDict["bgndFilePath"])
    #3. Tree name
    argDict["bgndTreeName"]=argDict["bgndTreeName"][0]    
    #TODO check if tree exists
    #4. energy branch name
    argDict["bgndEnergyBranchName"]=argDict["bgndEnergyBranchName"][0]    
    #TODO check if branch exists
    #5. Normalization option
    argDict["bgndNormalization"]=argDict["bgndNormalization"][0]    
    checkIfInputValid("bgndNormalization",argDict["bgndNormalization"],1,["fixed","branch"])
    #6. bgndNormalizationParameter
    if not argDict["bgndNormalization"]=="fixed":
        argDict["bgndNormalizationParameter"]=argDict["bgndNormalizationParameter"][0]
        #TODO check if branch exists
    else:
        argDict["bgndNormalizationParameter"]=float(argDict["bgndNormalizationParameter"][0])
    #7. bgndAdditionalCuts
    if not "bgndAdditionalCuts" in argDict:   
        argDict["bgndAdditionalCuts"]=[""]
    argDict["bgndAdditionalCuts"]=argDict["bgndAdditionalCuts"][0]
    
    #Check Source Data
    nSources=len(argDict["sourceNames"])
    nSims=len(argDict["simNames"])
    if not nSources==nSims:
        print("Need to supply same # of sim and source files!")
        print("Sources: "+str(nSources))
        print("Sims: "+str(nSims))
        sys.exit(-1)
    #1. sourceNames 
    #2. sourcePaths   
    for srcPath in argDict["sourcePaths"]:
        checkIfFileExists(srcPath)
    #3. sourceTreeName
    nTrees=len(argDict["sourceTreeName"])
    if not nTrees==nSources:
        print("Not enough source tree names provided, assuming same name for all sources")
        argDict["sourceTreeName"]=[argDict["sourceTreeName"][0] for i in range(0,nSources)]
    #4. sourceEnergyBranchName
    nEnergyBranches=len(argDict["sourceEnergyBranchName"])
    if not nEnergyBranches==nSources:
        print("Not enough source Energy Branches provided, assuming same for all sources")
        argDict["sourceEnergyBranchName"]=[argDict["sourceEnergyBranchName"][0] for i in range(0,nSources)]
    #5. sourceNormalization
    nNormalizations=len(argDict["sourceNormalization"])
    if not nNormalizations==nSources:
        print("Not enough source normalizations provided, assuming same for all sources")
        argDict["sourceNormalization"]=[argDict["sourceNormalization"][0] for i in range(0,nSources)]
    for i in range(0,nSources):
        checkIfInputValid("sourceNormalization",argDict["sourceNormalization"][i],1,["fixed","branch"])
    #6. sourceNormalizationParameter
    nNormalizationParameters=len(argDict["sourceNormalizationParameter"])
    if not nNormalizationParameters==nSources:
        print("Not enough source normalization parameters provided, assuming same for all sources")
        argDict["sourceNormalizationParameter"]=[argDict["sourceNormalizationParameter"][0] for i in range(0,nSources)]
    argDict["sourceNormalizationParameter"]=[float(argDict["sourceNormalizationParameter"][i]) if argDict["sourceNormalization"][i]=="fixed" else argDict["sourceNormalizationParameter"][i] for i in range(0,nSources)]
    #7. Source additional cuts
    if not "sourceAdditionalCuts" in argDict:
        argDict["sourceAdditionalCuts"]=["" for i in range(0,nSources)]
    nAdditionalCuts=len(argDict["sourceAdditionalCuts"])
    if not nAdditionalCuts==nSources:
        print("Not enough Additional Cuts parameters provided, assuming same for all sources")                  
        argDict["sourceAdditionalCuts"]=[argDict["sourceAdditionalCuts"][0] for i in range(0,nSources)]
        
    
    #Check Sim Data
    #1. simNames
    #2. simPaths   
    for simPath in argDict["simPaths"]:
        checkIfFileExists(simPath)
    #3. simTreeName
    nTrees=len(argDict["simTreeName"])
    if not nTrees==nSims:
        print("Not enough sim tree names provided, assuming same name for all sims")
        argDict["simTreeName"]=[argDict["simTreeName"][0] for i in range(0,nSims)]
    #4. sourceEnergyBranchName
    simBranchName=len(argDict["simBranchName"])
    if not nEnergyBranches==nSims:
        print("Not enough sim Energy Branches provided, assuming same for all sims")
        argDict["simBranchName"]=[argDict["simBranchName"][0] for i in range(0,nSims)]
    #5. simBranchUnits
    if not "simBranchUnits" in argDict:
        argDict["simBranchUnits"]=["MeV" for i in range(0,nSims)]
    nSimUnits=len(argDict["simBranchUnits"])
    if not nSimUnits==nSims:
        print("Not enough sim unit parameters provided, assuming same for all sims")                  
        argDict["simBranchUnits"]=[argDict["simBranchUnits"][0] for i in range(0,nSims)]
    for i in range(0,nSims):
        checkIfInputValid("simBranchUnits",argDict["simBranchUnits"][i],1,["keV","MeV"])
    #6. Sim additional cuts
    if not "simAdditionalCuts" in argDict:
        argDict["simAdditionalCuts"]=["" for i in range(0,nSims)]
    nAdditionalCuts=len(argDict["simAdditionalCuts"])
    if not nAdditionalCuts==nSims:
        print("Not enough Additional Cuts parameters provided, assuming same for all sims")                  
        argDict["simAdditionalCuts"]=[argDict["simAdditionalCuts"][0] for i in range(0,nSims)]
    
    return argDict  
    
#Does the parsing
def parseFile(fname):
    #Store input in a dictionary
    argDict={}
    argNamesList=[]
    argsList=[]
    inpFile=open(fname,"r")
    for lineNum,line in enumerate(inpFile):
        line=line.strip()
        #Ignore blank lines, comments
        if not line=="" and not line.startswith("#"):
            argName,args=parseLine(line)
            argDict[argName]=args
    argDict=checkInputForErrors(argDict)
    return argDict
    
#Custom class for holding data objects. Calling it calculates the normalization (if not provided)
class dataObject:
    def __init__(self,name,fileName,treeName,energyBranchName,normalizationOption,normalizationPar,additionalCuts,rangeMin,rangeMax):
        self.name=name
        print("Loading "+name+" data object")

        #Calculate normalization when we initialize, & load data into a numpy array
        df=ROOT.RDataFrame(treeName,fileName)
        if normalizationOption=="branch":
            npy=df.Filter(additionalCuts).AsNumpy(columns=[energyBranchName,normalizationPar])
            self.normalization=npy[normalizationPar][-1]-npy[normalizationPar][0]
        else:
            self.normalization=normalizationPar
            npy=df.Filter(additionalCuts).AsNumpy(columns=[energyBranchName])
        
        #Calculate # of counts in fit range, for generating expected # of counts
        rangeCut=energyBranchName+">="+str(rangeMin)+"&&"+energyBranchName+"<"+str(rangeMax)
        cut=rangeCut if additionalCuts=="" else rangeCut+"&&"+additionalCuts
        self.countsInFitRange = len(df.Filter(cut).AsNumpy(columns=[energyBranchName])[energyBranchName])
        print("Counts in fit range:"+str(self.countsInFitRange))
        self.data=npy[energyBranchName]
        print("Normalization = "+str(self.normalization))
    
#Custom class for holding sim objects. Converts energy to MeV if in keV
class simObject:
    def __init__(self,name,fileName,treeName,energyBranchName,branchUnits,maxSimEntries,additionalCuts):
        self.name=name
        print("Loading "+name+" sim object")
        #Calculate normalization when we initialize, & load data into a numpy array
        df=ROOT.RDataFrame(treeName,fileName)
        if branchUnits=="keV":
            npy=df.Filter(additionalCuts).AsNumpy(columns=[energyBranchName*1000.])
        else:
            npy=df.Filter(additionalCuts).AsNumpy(columns=[energyBranchName])
        if maxSimEntries>0:
            if maxSimEntries>len(npy[energyBranchName]):
                maxSimEntries=len(npy[energyBranchName])
        else:
            maxSimEntries=len(npy[energyBranchName])
                
        self.data=npy[energyBranchName][0:maxSimEntries]
        print("Loaded " +str(maxSimEntries)+" sim entries")



#Lambda fn for doing smearing
f = lambda x,alpha,beta,gamma: np.sqrt(np.power(alpha*x,2) + np.power(beta,2)*x + np.power(gamma,2))

#Input is a numpy array, output is a smeared array
def smearAndShiftArray(simArray,alpha,beta,gamma,slope,offset):
    #Smear each pt 20x
    valsToGen=20
    #Smear the MeV simulated data
    sigmas=f(simArray,alpha,beta,gamma)
    #Scale the simulation to match the data
    smearedArray=[]
    for i in range(0,valsToGen):
        smearedArray.append(slope*(np.random.normal(simArray,sigmas)-offset))

    #Make it a flat list
    return np.array(smearedArray).flatten()

def fillHist(hist,data):
    hist.Reset()
    w=np.full(data.size,1.)
    hist.FillN(data.size,data,w)
    return hist

def guessParams(simData,sourceHist):
    c1=ROOT.TCanvas("c1","c1")
    sourceHist.SetStats(0)
    sourceHist.SetMarkerColor(ROOT.kBlack)
    sourceHist.SetMarkerStyle(20)
    sourceHist.SetMarkerSize(3)
    
    sourceIntegral=sourceHist.Integral()
    
    simHist=sourceHist.Clone("simHist")
    simHist.SetLineColor(ROOT.kBlue)
    simHist.SetLineWidth(2)
    
    a=float(input("alpha="))
    b=float(input("beta="))
    g=float(input("gamma="))
    slope=float(input("slope="))
    offset=float(input("offset="))

    while(True):
        newData=smearAndShiftArray(simData,a,b,g,slope,offset)
        simHist=fillHist(simHist,newData)
        if simHist.Integral()==0:
            print("Something went wrong, your data is off the chart. Have the right units?")
        else:
            simHist.Scale(sourceIntegral/simHist.Integral())
        sourceHist.Draw("p")
        simHist.Draw("hist same")
        c1.Modified()
        c1.Update()
        par=input("Einter parameter to change [(a)lpha,(b)eta,(g)amma,(s)lope,(o)ffset,e(x)it)]:")
        if par=="x":
            sys.exit(-1)
        value=float(input("Input value:"))
        if par=="a": 
            a=value
        elif par=="b":
            b=value
        elif par=="g":
            g=value
        elif par=="s":
            slope=value
        elif par=="o":
            offset=value
    sys.exit(-1)
    
def plotResults(sourceHist,sourceCounts,bgndHist,bgndCounts,simHist,simCounts,modelHist,name):
    c1=ROOT.TCanvas("c1","c1")
    sourceHist.Scale(sourceCounts/sourceHist.Integral())
    sourceHist.SetStats(0)
    sourceHist.SetMarkerStyle(20)
    sourceHist.SetMarkerSize(1)
    sourceHist.SetLineColor(ROOT.kBlack)
    sourceHist.Draw("PE")
    
    bgndHist.Scale(bgndCounts/bgndHist.Integral())
    bgndHist.SetLineStyle(ROOT.kDashed)
    bgndHist.SetLineColor(ROOT.kBlack)
    bgndHist.SetLineWidth(2)
    bgndHist.Draw("Hist same")
    
    simHist.Scale(simCounts/simHist.Integral())
    simHist.SetLineColor(ROOT.kBlue)
    simHist.SetLineWidth(2)
    simHist.Draw("Hist same")
    
    modelHist.Reset()
    modelHist.Add(bgndHist,1)
    modelHist.Add(simHist,1)
    modelHist.SetLineColor(ROOT.kRed)
    modelHist.SetLineWidth(2)
    modelHist.Draw("Hist same")
    
    c1.Modified()
    c1.Update()
    c1.SaveAs(name+".pdf")
    
    
        
