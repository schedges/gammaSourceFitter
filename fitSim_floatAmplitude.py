# This version of the fitting code loops over energy resolution parameters,
# and energy calibration parameters, and fits the amplitude of a source
# pdf. The NLL value from the fit, along with the best fit parameters, are
# stored in a TTree. Profiles of the various parameters are also written
# to the output file.
import ROOT
import array
import numpy
import math
import gc
import os
import time
import multiprocessing

################
##Fit Settings##
################
# Either a "binned" background PDF or a "keys" background PDF.
bgndPdfType="binned"

# Number of bins to use for the simulation PDF, and for the bgnd PDF if
# above variable set to "binned.
# Note: if binning is too fine, you may have zero bins in the final PDF which
# can cause problems.
nBins=100

# Set to 1 to limit source entries to a maxSourceEntries, reduce to speed up
limitMaxSourceEntries=0
maxSourceEntries=150000

# Set to 1 to limit simulation entries to maxSimEntries, reduce to speed up
limitMaxSimEntries=1
maxSimEntries=150000

# How many smeared points to generate for each real point, makes a smoother PDF
# and helps to fit resolution parameters, but can slow things down
valsToGen=15

#####################
##Data I/O Settings##
#####################
# Pass in path to background file and source file, assuming they have the same
# tree structure
bgndFilename = "bgnd.root"
sourceFilename = "na22.root"

# Name of the data tree, branch names and types
# 'd'=double, "i"=int
dataTreeName = "leadTree"
# Only used for normalization purposes, can remove and normalize
# source and background data on your own
dataTreeTimeBranchName = "waveformStartTimeSecondsPart"
dataTreeTimeBranchType = 'd'
dataTreeChannelBranchName="channelID"
dataTreeChannelBranchType='i'
dataTreeEnergyBranchName ="integral"
dataTreeEnergyBranchType = 'd' #d=double, i=long

# Specify the channel number in data trees
dataTreeChannelNum=1

####################
##Sim I/O Settings##
#####################
# Pass in path to simulation file
simFilename = "na22sim.root"

# Specify name of the simulation tree
simTreeName = "combinedEnergyTree"
simTreeEnergyBranchName = "energy"
simTreeEnergyBranchType = 'd'
simTreeChannelBranchName = "cellNum"
simTreeChannelBranchType = 'd'

# Specify the channel number in the simulation tree
simTreeChannelNum=1

###########################
##Specify output filename##
###########################
outputFilename = "cal_output.root"

############################
##Parameter/Range Settings##
############################
####################################
# Range to fit in data tree energy units, try to avoid threshold regions
# or regions with sparse data if binning is fine
fitRangeMin=8000
fitRangeMax=120000
# Full range of data to import, should at least be the size of fit range
importLowerEnergyBound = 0
importUpperEnergyBound = 150000

#############################
##keV->ADC Start Parameters##
#############################
#Typical parameterization of light output as a fn of energy, L(E), >~40 keV:
#L = c * (E-E_0)
#where c=slope and E_0=offset, typically between 0 and 25 keV

#Specify slope min, slope max, and number of steps
slopeMin=50
slopeMax=70
numSlopeSteps=2

#Specify offset min, offset max, and number of steps
offsetMin=-10
offsetMax=25
numOffsetSteps=1

#########################
##Resolution parameters##
#########################
# For liquid scintillators, the typical parameterization of the energy resolution is as follows:
# FWHM/E = sqrt(alpha^2 + beta^2/E + gamma^2/E^2)
# Typically parameters vary with geometric size. Some example parameters can be found
# in the following papers:
# https://www.sciencedirect.com/science/article/pii/S0168900213003203
# https://www.sciencedirect.com/science/article/pii/0168900294906254
# https://www.sciencedirect.com/science/article/pii/016890029190331J
# https://www.sciencedirect.com/science/article/pii/S0168900217304795
# NOTE: The way I have things set up I'm pulling a factor of 1/sqrt(2.355) into
# the definition of this parameter. It seems like there are some variations in
# literature as to whether that 2.355 is there or not, but it should only matter
# when comparing your calibration parameters to other's.
alphaMin = 0.04
alphaMax = 0.12
numAlphaSteps = 2

betaMin = 0.05
betaMax = 3
numBetaSteps = 2

gammaMin = 0.00
gammaMax = 3
numGammaSteps=1

# For plastic scintillators, simpler form sometimes used FWHM/E = sqrt(a/E) =>
# sigma = 1/2.355 * sqrt(a*E)
# e.g. https://arxiv.org/abs/1004.3779 where a varied from 0.014 to 0.032 for
# those blocks. May need larger range for different sized detectors.
#
# Another group uses (sigma/E)^2 = A + B/E =>
# sigma = sqrt(A*E^2+B*E)
# and found A = 0.0162, B=0.0105 for their
# detectors: https://www.sciencedirect.com/science/article/pii/0029554X71903703
# They justified that the non-linear response for their plastic was zero, so there was no third term
#
# Note: You'll have to change a lot in the code if you change the resolution function

#Assume your PDFs will be within 20% of the expected number of counts
ampMin=0.8
ampMax=1.2

#Change if you are changing your calibration of resolution parameters. Last parameter should be "nll"
fitVarNames=["alpha","beta","gamma","slope","offset"]

################################################################################
##Assuming you followed the specified folder structure, you shouldn't need to###
##modify anything below here...#################################################
################################################################################
#* other than changing keyspdf setting or interpolation, or normalization



###################
##RooFit Settings##
###################
# Reduces warning messages
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

# Some fits failing with default (20), so I increased to 30
ROOT.RooAbsReal.defaultIntegratorConfig().getConfigSection("RooIntegrator1D").setRealValue("maxSteps",30);

# Set integration step sizes, useful for balancing speed vs. accuracy.
ROOT.RooAbsReal.defaultIntegratorConfig().setEpsAbs(1e-7)
ROOT.RooAbsReal.defaultIntegratorConfig().setEpsRel(1e-7)

################################
##Set up observable, fit range##
################################
energyVar = ROOT.RooRealVar("energyVar","energyVar",importLowerEnergyBound,importUpperEnergyBound)
energyVar.setRange("fitRange",fitRangeMin,fitRangeMax)
binning = ROOT.RooBinning(nBins,importLowerEnergyBound,importUpperEnergyBound,"binning")
energyVar.setBinning(binning)
# We'll apply this cut to the source and bgnd data sets solely for determining the expected
#counts in those PDFs for use in the fit
cut="energyVar>="+str(fitRangeMin)+"&&energyVar<="+str(fitRangeMax)
# Used a bunch of times, defining it once saves memory
argSet=ROOT.RooArgSet(energyVar)
argList=ROOT.RooArgList(energyVar)

############################
##Load background data set##
############################
bgndFile = ROOT.TFile(bgndFilename,"READ")
bgndTree = bgndFile.Get(dataTreeName)
# Make variables to read in branch quantities
dataTreeEnergy = array.array(dataTreeEnergyBranchType,[0])
dataTreeChannel = array.array(dataTreeChannelBranchType,[0])
dataTreeTime = array.array(dataTreeTimeBranchType,[0])
# Set branches
bgndTree.SetBranchAddress(dataTreeEnergyBranchName,dataTreeEnergy)
bgndTree.SetBranchAddress(dataTreeTimeBranchName,dataTreeTime)
bgndTree.SetBranchAddress(dataTreeChannelBranchName,dataTreeChannel)

# Create new data set to hold bgnd data
bgndDataSet = ROOT.RooDataSet("bgndDataSet","bgndDataSet",argSet)

# Get number of entries
nBgndEntries=bgndTree.GetEntries()
print("Found "+str(nBgndEntries)+" in bgnd tree")

# Step through bgnd tree, adding data if meets energy and channel criteria
for entry in range(0,nBgndEntries):
  bgndTree.GetEntry(entry)
  # Check if this belongs to the channel number we're fitting data for
  if dataTreeChannel[0] == dataTreeChannelNum:
    # Check if the value is in the import range
    if (dataTreeEnergy[0] >= importLowerEnergyBound) and (dataTreeEnergy[0] <= importUpperEnergyBound):
      # Set value of observable
      energyVar.setVal(dataTreeEnergy[0])
      # Add to data set
      bgndDataSet.add(argSet)

#Get length or run for normalization
#Get first entry
bgndTree.GetEntry(0)
startTime = dataTreeTime[0]
#Get last entry
bgndTree.GetEntry(nBgndEntries-1)
endTime = dataTreeTime[0]
#Calculate run length
bgndLengthOfRun = endTime-startTime

#Get number of counts in fit range
reducedDataSet = bgndDataSet.reduce(ROOT.RooFit.Cut(cut))
bgndCountsInFitRange = reducedDataSet.numEntries()
print("Found "+str(bgndCountsInFitRange)+" bgnd entries in the fit range")

#Make background pdf
if bgndPdfType=="binned":
  print("Making binned bgnd Pdf\n")
  (bgndDataSet.get().find("energyVar")).setBins(nBins)
  bgndDataHist = bgndDataSet.binnedClone()
  bgndDataPdf = ROOT.RooHistPdf("bgndDataPdf","bgndDataPdf",argSet,bgndDataHist,1) #1 specifies interpolation order
else:
  print("Making bgnd RooKeysPdf\n")
  bgndDataPdf = ROOT.RooKeysPdf("bgndDataPdf","bgndDataPdf",energyVar,bgndDataSet,ROOT.RooKeysPdf.NoMirror,1.2)

########################
##Load source data set##
########################
sourceFile = ROOT.TFile(sourceFilename,"READ")
sourceTree = sourceFile.Get(dataTreeName)
#Set branches
sourceTree.SetBranchAddress(dataTreeEnergyBranchName,dataTreeEnergy)
sourceTree.SetBranchAddress(dataTreeTimeBranchName,dataTreeTime)
sourceTree.SetBranchAddress(dataTreeChannelBranchName,dataTreeChannel)

#Create RooDataSet for holding source data
sourceDataSet = ROOT.RooDataSet("sourceDataSet","sourceDataSet",argSet)

#Get number of entries
nSourceEntries=sourceTree.GetEntries()
print("Found "+str(nSourceEntries)+" entries in source data set")

if limitMaxSourceEntries==1:
  if nSourceEntries>maxSourceEntries:
    nSourceEntries=maxSourceEntries
    print("Limiting source to "+str(nSourceEntries))

#Step through source data set, adding entries if meet the appropriate energy and channel criteria
for entry in range(0,nSourceEntries):
  sourceTree.GetEntry(entry)
  if dataTreeChannel[0] == dataTreeChannelNum:
    if (dataTreeEnergy[0] >= importLowerEnergyBound) and (dataTreeEnergy[0] <= importUpperEnergyBound):
      energyVar.setVal(dataTreeEnergy[0])
      sourceDataSet.add(argSet)

#Get length or run for normalization
#Get first entry
sourceTree.GetEntry(0)
startTime = dataTreeTime[0]
#Get last entry
sourceTree.GetEntry(nSourceEntries-1)
endTime = dataTreeTime[0]
#Calculate length of run
sourceLengthOfRun = endTime-startTime

#Get number of counts in fit range
reducedDataSet = sourceDataSet.reduce(ROOT.RooFit.Cut(cut))
sourceCountsInFitRange = reducedDataSet.numEntries()
print("Found "+str(sourceCountsInFitRange)+" source entries in the fit range\n")

#Make binned source PDF
(sourceDataSet.get().find("energyVar")).setBins(nBins)
sourceDataHist = sourceDataSet.binnedClone()
srcDataHistPdf = ROOT.RooHistPdf("srcDataHistPdf","srcDataHistPdf",argSet,sourceDataHist,0) #1 specifies interpolation order

#################################################
##Normalize bgnd counts to length of source run##
#################################################
#Scale length of background run to match source run
scaleFactor = float(sourceLengthOfRun)/float(bgndLengthOfRun)
scaledBgndEntries = bgndCountsInFitRange*scaleFactor
#Make a RooRealVar to hold this, and set expected bgnd counts constant
scaledBgndEntriesVar = ROOT.RooRealVar("scaledBgndEntriesVar","scaledBgndEntriesVar",scaledBgndEntries)
scaledBgndEntriesVar.setConstant(1)
print("Length of bgnd run: "+str(bgndLengthOfRun))
print("Length of source run: "+str(sourceLengthOfRun))
print("Scale factor: "+str(scaleFactor))
print("Scaled bgnd counts in fit range: "+str(scaledBgndEntries))
print("Expected source counts in fit range: "+str(sourceCountsInFitRange-scaledBgndEntries)+"\n")
expectedSourceCounts=sourceCountsInFitRange-scaledBgndEntries

#####################
##Load Sim Data Set##
#####################
#Open sim file
simFile = ROOT.TFile(simFilename,"READ")
#Get simTree from file
simTree = simFile.Get(simTreeName)
#Make variables to hold sim tree quantities
simTreeEnergy = array.array(simTreeEnergyBranchType,[0])
simTreeChannel = array.array(simTreeChannelBranchType,[0])
#Load branches
simTree.SetBranchAddress(simTreeEnergyBranchName,simTreeEnergy)
simTree.SetBranchAddress(simTreeChannelBranchName,simTreeChannel)

#Get number of entries
nSimEntries=simTree.GetEntries()
print("Found "+str(nSimEntries)+" entries in simulation tree")
if limitMaxSimEntries==1:
  if nSimEntries>maxSimEntries:
    nSimEntries=maxSimEntries
    print("Limiting sim to "+str(nSimEntries))

#Store sim data in a numpy array for fast smearing
simList=[]
for entry in range(0,nSimEntries):
  simTree.GetEntry(entry)
  simList.append(simTreeEnergy[0])
simArray=numpy.array(simList)

#####################################
##Make tree for storing fit results##
#####################################
outputFile = ROOT.TFile(outputFilename,"RECREATE")
fitTree = ROOT.TTree("fitTree","Results from fit")

alpha=array.array('d',[0])
beta=array.array('d',[0])
gamma=array.array('d',[0])
slope=array.array('d',[0])
offset=array.array('d',[0])
fitNll=array.array('d',[0])
bgndCounts=array.array('d',[0])
bgndCounts[0]=scaledBgndEntries
sourceCounts=array.array('d',[0])
sourceCountsError=array.array('d',[0])
fitStatus=array.array('i',[0])

fitTree.Branch('alpha',alpha,'alpha/D')
fitTree.Branch('beta',beta,'beta/D')
fitTree.Branch('gamma',gamma,'gamma/D')
fitTree.Branch('slope',slope,'slope/D')
fitTree.Branch('offset',offset,'offset/D')
fitTree.Branch('fitNll',fitNll,'fitNll/D')
fitTree.Branch('bgndCounts',bgndCounts,'bgndCounts/D')
fitTree.Branch('sourceCounts',sourceCounts,'sourceCounts/D')
fitTree.Branch('sourceCountsError',sourceCountsError,'sourceCountsError/D')
fitTree.Branch('fitStatus',fitStatus,'fitStatus/I')

################################
##Making parameter sets to fit##
################################
#Make a set of all possible combinations of slope, offset, and resolution parameters
alphaSpace=numpy.linspace(alphaMin,alphaMax,numAlphaSteps)
betaSpace=numpy.linspace(betaMin,betaMax,numBetaSteps)
gammaSpace=numpy.linspace(gammaMin,gammaMax,numGammaSteps)
slopeSpace=numpy.linspace(slopeMin,slopeMax,numSlopeSteps)
offsetSpace=numpy.linspace(offsetMin,offsetMax,numOffsetSteps)
 
#Combine to make a parameter set we'll iterate through
parametersToFit=[]
for i in alphaSpace:
  for j in betaSpace:
    for k in gammaSpace:
      for l in slopeSpace:
        for m in offsetSpace:
          parametersToFit.append([i,j,k,l,m])

#Efficient way to apply resolution to an array of values
f = lambda x,alpha,beta,gamma: numpy.sqrt( numpy.power(alpha*x,2) + numpy.power(beta,2)*x + numpy.power(gamma,2) )

#Set to 1 to plot at the end
PLOT=0

###########
##Fitting##
###########
def smearAndFit(parSet):
  global bestParSet

  #Set Parameters
  alpha=parSet[0]
  beta=parSet[1]
  gamma=parSet[2]
  slope=parSet[3]
  offset=parSet[4]

  #Reset source amplitude var to expected counts
  sourceCountsVar = ROOT.RooRealVar("sourceCountsVar","sourceCountsVar",expectedSourceCounts,ampMin*expectedSourceCounts,ampMax*expectedSourceCounts)
  sourceCountsVar.setVal(expectedSourceCounts)

  #Make array of sigmas from res params and sim values
  sigmas=f(simArray,alpha,beta,gamma)

  #Generate valsToGen random values for every entry in sim
  smearedArrayList=[]
  for i in range(0,valsToGen):
    smearedArrayList.append(slope*(numpy.random.normal(simArray,sigmas)-offset))

  #Make it a flat list
  smearedArrayArray=numpy.array(smearedArrayList)
  flatArray=smearedArrayArray.flatten()

  #Numpy array->TH1->RooDataHist->RooHistPdf
  #~0.03 seconds!!!
  smearedSimHist=ROOT.TH1D("smearedSimHist","smearedSimHist",nBins,importLowerEnergyBound,importUpperEnergyBound)
  weights=numpy.full(flatArray.size,1.)
  smearedSimHist.FillN(flatArray.size,flatArray,weights)
  smearedSimDataHist=ROOT.RooDataHist("smearedSimDataHist","smearedSimDataHist",argList,smearedSimHist)
  simPdf = ROOT.RooHistPdf("simPdf","simPdf",argSet,smearedSimDataHist,0) #1 specifies interpolation order

  ##Make Model
  pdfList = ROOT.RooArgList(bgndDataPdf,simPdf)
  ampVarList = ROOT.RooArgList(scaledBgndEntriesVar,sourceCountsVar)
  model = ROOT.RooAddPdf("model","model",pdfList,ampVarList)
  model.fixCoefRange("fitRange")

  #Fit
  res = model.fitTo(sourceDataSet,
    ROOT.RooFit.Extended(1),
    ROOT.RooFit.Range("fitRange"),
    ROOT.RooFit.SumCoefRange("fitRange"),
    ROOT.RooFit.Verbose(0),
    ROOT.RooFit.PrintLevel(-1),
    ROOT.RooFit.Save(True)
  )
  
  #Set rest of quantities in ttree and fill
  fitNll = res.minNll();
  fitStatus = res.status()
  sourceCounts = sourceCountsVar.getVal()
  sourceCountsError = sourceCountsVar.getError()
  
  if PLOT==1:
    try:
      c1
    except NameError:
      c1=ROOT.TCanvas("c1","c1")
      
    #Reduce integrator for plotting, massively speeds things up
    ROOT.RooAbsReal.defaultIntegratorConfig().setEpsAbs(1e-5)
    ROOT.RooAbsReal.defaultIntegratorConfig().setEpsRel(1e-5)
     
    frame = energyVar.frame(importLowerEnergyBound,importUpperEnergyBound,nBins)
    frame.SetTitle("Alpha="+str(alpha)+", Beta="+str(beta)+", Gamma="+str(gamma)+", Slope="+str(slope)+", Offset="+str(offset))

    #Plot source data
    sourceDataSet.plotOn(
      frame,
      ROOT.RooFit.Name("Source"),
      ROOT.RooFit.MarkerColor(1),
      ROOT.RooFit.FillColor(0),
    )
    #Plot components
    model.plotOn(
      frame,
      ROOT.RooFit.Name("Bgnd"),
      ROOT.RooFit.Components("bgndDataPdf"),
      ROOT.RooFit.LineColor(ROOT.kSolid),
      ROOT.RooFit.FillColor(0),
      ROOT.RooFit.ProjWData(sourceDataSet)
    )
    model.plotOn(
      frame,ROOT.RooFit.Name("Sim"),
      ROOT.RooFit.Components("simPdf"),
      ROOT.RooFit.LineColor(ROOT.kRed),
      ROOT.RooFit.FillColor(0),
      ROOT.RooFit.ProjWData(sourceDataSet),
      ROOT.RooFit.AddTo("Bgnd")
    )

    #Draw
    frame.Draw()

    #Add legend
    leg = ROOT.TLegend(0.65,0.65,0.95,0.92);
    sourceObj = frame.findObject("Source");
    bgndObj = frame.findObject("Bgnd");
    simObj = frame.findObject("Sim");
    leg.AddEntry(sourceObj,"Source","P")
    leg.AddEntry(bgndObj,"Background","L")
    leg.AddEntry(simObj,"Sim","L")
    leg.Draw("same")
  
    #Draw
    c1.SetLogy()
    c1.Modified()
    c1.Update()
    c1.SaveAs("bestFit_ch"+str(dataTreeChannelNum)+".pdf")
  
  #Memory management
  sourceCountsVar.Delete()
  del sourceCountsVar
  pdfList.Delete()
  del pdfList
  ampVarList.Delete()
  del ampVarList
  smearedSimHist.Delete()
  del smearedSimHist
  smearedSimDataHist.Delete()
  del smearedSimDataHist
  simPdf.Delete()
  del simPdf
  #model.Delete()
  #del model
  #res.Delete()
  #del res
  gc.collect()
  
  return alpha,beta,gamma,slope,offset,fitNll,fitStatus,sourceCounts,sourceCountsError
 
##############
##Do the fits#
##############
#We'll normalize the NLLs so the lowest is zero. To do this we search for the lowest
#value and store in bestNll
bestNll=numpy.inf
bestTheta=""
results=[]
averageTime=0
print("")
stepNum=0
for parSet in parametersToFit:
  startTime=time.time()
  theta=smearAndFit(parSet)
  results.append(theta)
  endTime=time.time()-startTime
  
  stepNum+=1
  averageTime = (averageTime*(stepNum-1)+endTime)*1./stepNum
  if stepNum!=len(parametersToFit):
    print("Step {} of {}, elapsed time/estimated time: {:.1f}/{:.1f} sec., average iteration time : {:0.1f} sec.".format(stepNum,len(parametersToFit),averageTime*stepNum,len(parametersToFit)*averageTime,averageTime),end="\r")
  else:
    print("\nFinished fitting in "+str(averageTime*stepNum)+" seconds!")
    

  currentNLL=theta[len(fitVarNames)]
  if currentNLL<bestNll:
    bestNll=currentNLL
    bestTheta=theta

print("Best nll: "+str(bestNll)+" with par set: "+str(bestTheta))

#Fill TTree
for theta in results:
  alpha[0]=theta[0]
  beta[0]=theta[1]
  gamma[0]=theta[2]
  slope[0]=theta[3]
  offset[0]=theta[4]
  fitNll[0]=theta[5]-bestNll
  fitStatus[0]=theta[6]
  sourceCounts[0]=theta[7]
  sourceCountsError[0]=theta[8]
  
  fitTree.Fill()

#####################
##Gen profile plots##
#####################
shortenedResults=[]
for result in results:
  element=result[0:len(fitVarNames)+1]
  shortenedResults.append(element)
  
graphs=[]
#Step through every component but NLL
for parNum in range(0, len(fitVarNames)):

	graphs.append(ROOT.TGraph())
	graphs[-1].SetTitle(fitVarNames[parNum])
	graphs[-1].SetName(fitVarNames[parNum])
	graphs[-1].GetXaxis().SetTitle(fitVarNames[parNum])
	graphs[-1].GetYaxis().SetTitle("NLL")
	nPoints=0
	
	#Make list of unique values of this parameter
	uniqueElementList=[]
	for parSet in results:
		if parSet[parNum] not in uniqueElementList:
			uniqueElementList.append(parSet[parNum])

	#Now calculate min nll for each value 
	minNlls=[]
	for uniqueElement in uniqueElementList:
		nlls=[]
		for parSet in results:
			if parSet[parNum]==uniqueElement:
				nlls.append(parSet[-1])

		minNll = min(nlls)
		minNlls.append(minNll)

	#Offset so lowerest NLL occurs at zero
	for uniqueElement in uniqueElementList:
		graphs[-1].SetPoint(nPoints,uniqueElement,minNlls[nPoints])
		nPoints+=1

	graphs[-1].SetMarkerStyle(20)

PLOT=1
smearAndFit(bestTheta)

#####################
#Write trees to file#
#####################
outputFile.cd()
fitTree.Write()
for graph in graphs:
	graph.Write()
outputFile.Close()
