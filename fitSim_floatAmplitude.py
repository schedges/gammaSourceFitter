#This version of the fitting code loops over values of different fitting parameters in discrete steps,
#fits the ampltitude of the simulation PDF, and returns the NLL value. This utilizes a binned simulation
#and bgnd PDFs. Care should be taken to test the effect of binning on parameter values.
#
#A separate version to sample the parameter space using Markov-Chain Monte Carlo is being developed.
#
#This code outputs three things to the specified output file
#	fitTree: containing fit parameters for each fit
#	bestFitTree: containing fit parameters for fit with lowest NLL only
#	canvases of fits: Currently only saves fits it finds with a lower NLL than the current lowest NLL value,
#		so last canvas saved should correspond to the best fit.
# profile likelihoods of each parameters--to be used in determining ranges, fit errors
#
#Requirements:
#	Source/bgnd data in root TTree with:
#		- Branch corresponding to energy
#		- Branch corresponding to time (otherwise you need to manually normalize data)
#		- Branch corresponding to channel number (can modify the code if only fitting one channel)
#	Simulation data TTree with:
#		- Branch corresponding to energy (keV)
#		- Branch corresponding to channel (can optionally modify code if only fitting one channel)
#
#Final notes:
#		- Fitting a spectrum with a single gamma line has not worked well for me. It
#		seems fairly common in the literature for sims to underpredict the low energy 
#		part of the spectrum. As a result, the fits tend to favor very large negative 
#		offsets to counteract this. This code could be modified to do a simultaneous 
#		fit with another source to keep the offset from growing unphysically large, a
#		more complex gamma source could be used, or external constraints could be applied
#		to the offset. 
import ROOT
import array
import numpy
import math
import gc

#Specify number of bins for PDFs, over the entire import range (importLowerEnergyBound, importUpperEnergyBound).
#NOTE: If this is too high, there may be zero bins in the bgnd pdf that then 
#cause problems when evaluating the NLL. 
#This is only used for plotting fits if using a bgnd and simulation Keys Pdf
nBins=100

#If you use these, you should plot the keys PDF to ensure it matches the data. You 
#may also need to manually edit the KeysPdf generation lines to change the smoothing
useBgndKeysPdf=1
useSimKeysPdf=1

#Specify how many threads to use for fit
numCPUForFits=1

#Decrease fraction to speed things up, useful for tests
fractionOfDataToUse=1.0

#How many smeared points to generate for each real point, makes a smoother PDF
#without requiring increase to number of sim events
valsToGen=1

#Canvas for plotting 
c1 = ROOT.TCanvas("c1","c1")



##################################
##Specify data format/names here##
##################################
#Specify the channel number in data trees
dataTreeChannelNum=0

#Specify name of bgnd file and source file (and path if not in current directory).
#This code assumes the bgnd data set and source data set have the same tree name, branches, etc.
bgndFilename = "bgnd.root"
sourceFilename = "na22.root"

#Specify name of data tree
dataTreeName = "BDTree"
#Specify name and type of branch corresponding to energy
dataTreeEnergyBranchName ="integral"
dataTreeEnergyBranchType = 'd' #d=double, i=long
dataTreeEnergy = array.array(dataTreeEnergyBranchType,[0])
#Specify name and type of branch corresponding to channelNum
dataTreeChannelBranchName="channelNum"
dataTreeChannelBranchType='i'
dataTreeChannel = array.array(dataTreeChannelBranchType,[0])
#Specify name and type of branch corresponding to run time (for normalization)
dataTreeTimeBranchName = "waveformTime"
dataTreeTimeBranchType = 'd'
dataTreeTime = array.array(dataTreeTimeBranchType,[0])



#################################
##Specify sim format/names here##
#################################
#For sim trees with multiple channels, specify the channel number here
simTreeChannelNum=0

#Specify name of simulation file (and path if not in current directory)
simFilename = "sim-na22-1e7.root"

#Specify name of the simulation tree
simTreeName = "BD_cal"
#Specify name of the energy branch in the simulation tree
simTreeEnergyBranchName = "energy"
simTreeEnergyBranchType = 'd'
simTreeEnergy = array.array(simTreeEnergyBranchType,[0])
#Specify the name of channel number branch in the simulation tree
simTreeChannelBranchName = "BD_id"
simTreeChannelBranchType = 'i'
simTreeChannel = array.array(simTreeChannelBranchType,[0])



##############
#Output name##
##############
outputFilename = "ch"+str(dataTreeChannelNum)+"_cal.root"



####################################
##Specify data import / fit ranges##
####################################
#Specify the range to fit over (in terms of dataTreeEnergy variable).
#NOTE: you probably want the minimum value to be above the threshold region, otherwise
#threshold effects can upset the fit
fitRangeMin=1250
fitRangeMax=25000
#Specify range for importing energy variables from both dataTrees and simTrees.
#NOTE: this should be larger than the fit range, as data from outside the fit range
#can be brought into it when applying smearing. I've also had problems when the bounds
#includes 0, so I start just above that.
#NOTE: These values are inclusive.
importLowerEnergyBound = 1
importUpperEnergyBound = 30000



#############################
##keV->ADC Start Parameters##
#############################
#Typical parameterization of light output as a fn of energy, L(E), above 40 keV, is
#L = c * (E-E_0)
#where c=slope, the conversion factor, and E_0=offset,  a small offset,
#usually between 0 and 25 keV

#Specify slope range and number of steps to include between that range here. The 
#units should be ADC/keV
slopeMin=14
slopeMax=24
numSlopeSteps=6

#Specify offset here. Note that because of the way I define offset later, I 
#think this should be a negative keV value, but I've had some fits favor positive
#offsets. The units of it should be keV.
offsetMin=-10
offsetMax=10
numOffsetSteps=5


#########################
##Resolution parameters##
#########################
#For liquid scintillators, the typical parameterization of the energy resolution is as follows:
#FWHM/E = sqrt(alpha^2 + beta^2/E + gamma^2/E^2)
#Typically parameters vary with geometric size. Some example parameters can be found
#in the following papers:
# https://www.sciencedirect.com/science/article/pii/S0168900213003203
# https://www.sciencedirect.com/science/article/pii/0168900294906254
# https://www.sciencedirect.com/science/article/pii/016890029190331J
# https://www.sciencedirect.com/science/article/pii/S0168900217304795
#NOTE: The way I have things set up I'm pulling a factor of 1/sqrt(2.355) into 
# the definition of this parameter. It seems like there are some variations in 
# literature as to whether that 2.355 is there or not, but it should only matter
# when comparing your calibration parameters to other's.
alphaMin = 0.04
alphaMax = 0.12
numAlphaSteps = 5

betaMin = 0.05
betaMax = 0.20
numBetaSteps = 6

gammaMin = 0.00
gammaMax = 0.06
numGammaSteps= 4

#For printing purposes only
totalSteps = numSlopeSteps*numOffsetSteps*numAlphaSteps*numBetaSteps*numGammaSteps


#For plastic scintillators, simpler form sometimes used FWHM/E = sqrt(a/E)
# => sigma = 1/2.355 * sqrt(a*E)
#e.g. https://arxiv.org/abs/1004.3779 where a varied from 0.014 to 0.032 for 
#those blocks. May need larger range for different sized detectors.
#
#Another group uses (sigma/E)^2 = A + B/E
#=> sigma = sqrt(A*E^2+B*E)
#and found A = 0.0162, B=0.0105 for their
#detectors: https://www.sciencedirect.com/science/article/pii/0029554X71903703
#They justified that the non-linear response for their plastic was zero, so there was no third term
#
#NOTE: You will need to change how the profiles are generated if you use a different energy 
#	resolution function


###################
##RooFit Settings##
###################
#Reduces warning messages. I suggest leaving these on at first to make sure there's no serious
#warning/error messages
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

#Controls the number of steps used for integration. I found fits sometimes failing
# with the default so I increased it to 30
ROOT.RooAbsReal.defaultIntegratorConfig().getConfigSection("RooIntegrator1D").setRealValue("maxSteps",30);

#Set integration step sizes, useful for balancing speed vs. accuracy. I believe 1e-7 is default
ROOT.RooAbsReal.defaultIntegratorConfig().setEpsAbs(1e-7)
ROOT.RooAbsReal.defaultIntegratorConfig().setEpsRel(1e-7)



########################MAIN CODE STARTS HERE###########################
#To get started, you shouldn't need to adjust anything below this line##



################################
##Set up observable, fit range##
################################
energyVar = ROOT.RooRealVar("energyVar","energyVar",importLowerEnergyBound,importUpperEnergyBound)
energyVar.setRange("fitRange",fitRangeMin,fitRangeMax)
#If either sim or source PDF is binned, define binning variable 
if useBgndKeysPdf==0 or useSimKeysPdf==0:
	binning = ROOT.RooBinning(nBins,importLowerEnergyBound,importUpperEnergyBound,"binning")
	energyVar.setBinning(binning)
#We'll apply this cut to the source and bgnd data sets solely for determining the expected
#counts in those PDFs for use in the fit
cut="energyVar>="+str(fitRangeMin)+"&&energyVar<="+str(fitRangeMax)
#Wrap energyVar into a RooArgSet, do it here so we don't create a bunch of these
#as it will eat up memory defining it inside RooFit commands 
argSet=ROOT.RooArgSet(energyVar)



############################
##Load background data set## 
############################
#Open bgnd file
bgndFile = ROOT.TFile(bgndFilename,"READ")
#Get tree from bgnd file
bgndTree = bgndFile.Get(dataTreeName)
#Set branches
bgndTree.SetBranchAddress(dataTreeEnergyBranchName,dataTreeEnergy)
bgndTree.SetBranchAddress(dataTreeTimeBranchName,dataTreeTime)
bgndTree.SetBranchAddress(dataTreeChannelBranchName,dataTreeChannel)

#Create new data set to hold bgnd data
bgndDataSet = ROOT.RooDataSet("bgndDataSet","bgndDataSet",ROOT.RooArgSet(energyVar))

#Get number of entries
nBgndEntries=int(bgndTree.GetEntries()*fractionOfDataToUse)
print("Found "+str(nBgndEntries)+" in bgnd tree")

#Step through bgnd tree, adding data if meets energy and channel criteria
for entry in range(0,nBgndEntries):
	bgndTree.GetEntry(entry)
	#Check if this belongs to the channel number we're fitting data for
	if dataTreeChannel[0] == dataTreeChannelNum:
		#Check if the value is in the import range
		if (dataTreeEnergy[0] >= importLowerEnergyBound) and (dataTreeEnergy[0] <= importUpperEnergyBound):
			#Set value of observable
			energyVar.setVal(dataTreeEnergy[0])
			#Add to data set
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
if useBgndKeysPdf==0:
	(bgndDataSet.get().find("energyVar")).setBins(nBins)
	bgndDataHist = bgndDataSet.binnedClone()
	bgndDataPdf = ROOT.RooHistPdf("bgndDataPdf","bgndDataPdf",argSet,bgndDataHist,1) #1 specifies interpolation order
else:
	bgndDataPdf = ROOT.RooKeysPdf("bgndDataPdf","bgndDataPdf",energyVar,bgndDataSet,ROOT.RooKeysPdf.NoMirror,0.7)



########################
##Load source data set## 
########################
#Open source file
sourceFile = ROOT.TFile(sourceFilename,"READ")
#Get source tree
sourceTree = sourceFile.Get(dataTreeName)
#Set branches
sourceTree.SetBranchAddress(dataTreeEnergyBranchName,dataTreeEnergy)
sourceTree.SetBranchAddress(dataTreeTimeBranchName,dataTreeTime)
sourceTree.SetBranchAddress(dataTreeChannelBranchName,dataTreeChannel)

#Create RooDataSet for holding source data
sourceDataSet = ROOT.RooDataSet("sourceDataSet","sourceDataSet",argSet)

#Get number of entries
nSourceEntries=int(sourceTree.GetEntries()*fractionOfDataToUse)
print("Found "+str(nSourceEntries)+" entries in source data set")

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

#Uses binned source data can speed up the fits, but definitely shouldn't be used
#if you're using Keys PDFs for bgnd and sim because you want an unbinned fit. 
#Specify whether you want to fit to the sourceDataSet or sourceDataHist in the fitting line
(sourceDataSet.get().find("energyVar")).setBins(nBins)
sourceDataHist = sourceDataSet.binnedClone()
srcDataHistPdf = ROOT.RooHistPdf("srcDataHistPdf","srcDataHistPdf",ROOT.RooArgSet(energyVar),sourceDataHist,1) #1 specifies interpolation order



#####################
##Load Sim Data Set##
#####################
#Open sim file
simFile = ROOT.TFile(simFilename,"READ")
#Get simTree from file
simTree = simFile.Get(simTreeName)
#Load branches
simTree.SetBranchAddress(simTreeEnergyBranchName,simTreeEnergy)
simTree.SetBranchAddress(simTreeChannelBranchName,simTreeChannel)

#Create data set for sim
simDataSet = ROOT.RooDataSet("simDataSet","simDataSet",argSet)

#Get number of entries
nSimEntries=simTree.GetEntries()
print("Found "+str(nSimEntries)+" entries in simulation tree\n")

#Step through, loading entries to data set if meet energy and channel criteria
for entry in range(0,nSimEntries):
	simTree.GetEntry(entry)
	if simTreeChannel[0]==simTreeChannelNum:
		if (simTreeEnergy[0] >= importLowerEnergyBound) and (simTreeEnergy[0] <= importUpperEnergyBound):
			energyVar.setVal(simTreeEnergy[0])
			simDataSet.add(argSet)

#Make data set for holding smeared data (smearing done in loop)
smearedSimDataSet = ROOT.RooDataSet("smearedSimDataSet","smearedSimDataSet",ROOT.RooArgSet(energyVar))



#################################################
##Normalize bgnd counts to length of source run##
#################################################
#Scale length of background run to match source run
scaleFactor = float(sourceLengthOfRun)/float(bgndLengthOfRun)
scaledBgndEntries = bgndCountsInFitRange*scaleFactor
#Make a RooRealVar to hold this, and set expected bgnd counts constant
scaledBgndEntriesVar = ROOT.RooRealVar("scaledBgndEntriesVar","scaledBgndEntriesVar",scaledBgndEntries)
scaledBgndEntriesVar.setConstant()
print("Length of bgnd run: "+str(bgndLengthOfRun))
print("Length of source run: "+str(sourceLengthOfRun))
print("Scale factor: "+str(scaleFactor))
print("Scaled bgnd counts in fit range: "+str(scaledBgndEntries))
print("Expected source counts in fit range: "+str(sourceCountsInFitRange-scaledBgndEntries))



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
fitNLL=array.array('d',[0])
fitBgndCounts=array.array('d',[0])
fitSourceCounts=array.array('d',[0])
fitSourceCountsError=array.array('d',[0])
fitStatus=array.array('i',[0])

fitTree.Branch('alpha',alpha,'alpha/D')
fitTree.Branch('beta',beta,'beta/D')
fitTree.Branch('gamma',gamma,'gamma/D')
fitTree.Branch('slope',slope,'slope/D')
fitTree.Branch('offset',offset,'offset/D')
fitTree.Branch('fitNLL',fitNLL,'fitNLL/D')
fitTree.Branch('fitBgndCounts',fitBgndCounts,'fitBgndCounts/D')
fitTree.Branch('fitSourceCounts',fitSourceCounts,'fitSourceCounts/D')
fitTree.Branch('fitSourceCountsError',fitSourceCountsError,'fitSourceCountsError/D')
fitTree.Branch('fitStatus',fitStatus,'fitStatus/I')



#############################################################
##Make a separate tree just for holding the best fit values##
#############################################################
bestFitTree = ROOT.TTree("bestFitTree","Best fit")
bestFitTree.Branch('alpha',alpha,'alpha/D')
bestFitTree.Branch('beta',beta,'beta/D')
bestFitTree.Branch('gamma',gamma,'gamma/D')
bestFitTree.Branch('slope',slope,'slope/D')
bestFitTree.Branch('offset',offset,'offset/D')
bestFitTree.Branch('fitNLL',fitNLL,'fitNLL/D')
bestFitTree.Branch('fitBgndCounts',fitBgndCounts,'fitBgndCounts/D')
bestFitTree.Branch('fitSourceCounts',fitSourceCounts,'fitSourceCounts/D')
bestFitTree.Branch('fitSourceCountsError',fitSourceCountsError,'fitSourceCountsError/D')
bestFitTree.Branch('fitStatus',fitStatus,'fitStatus/I')

bestNLL=10000000 #Set to something large so hopefully first fit has a lower nll and we plot it
bestSlope=0
bestOffset=0
bestAlpha=0
bestBeta=0
bestGamma=0
bestBgndCounts=0
bestSourceCounts=0
bestSourceCountsError=0
bestFitStatus=0


###########
##Fitting##
###########
stepNum=0
for testAlpha in numpy.linspace(alphaMin,alphaMax,numAlphaSteps):
	for testBeta in numpy.linspace(betaMin,betaMax,numBetaSteps):
		for testGamma in numpy.linspace(gammaMin,gammaMax,numGammaSteps):
			for testSlope in numpy.linspace(slopeMin,slopeMax,numSlopeSteps):
				for testOffset in numpy.linspace(offsetMin,offsetMax,numOffsetSteps):
			
					stepNum=stepNum+1
					print("\n\nOn step "+str(stepNum)+" of "+str(totalSteps))
					print("Alpha="+str(testAlpha)+", beta="+str(testBeta)+", gamma="+str(testGamma)+", slope="+str(testSlope)+", offset="+str(testOffset))
					
					
					
					#####################################
					##Get expected amplitude of sim PDF##
					#####################################
					expectedSourceCounts=sourceCountsInFitRange-scaledBgndEntries
					#Assume we're within 50% of the true values
					sourceCountsVar = ROOT.RooRealVar("sourceCountsVar","sourceCountsVar",expectedSourceCounts,0.5*expectedSourceCounts,1.5*expectedSourceCounts)
					#Ensure sourceCountsVar can float
					sourceCountsVar.setConstant(0)
					
					
					
					#####################
					##Smear simulations##
					#####################
					#Clear all exisitng entries in smeared sim data set
					smearedSimDataSet.reset()
					
					#Step through sim data set, smearing
					nSimEntries=simDataSet.numEntries()
					for entry in range(0,nSimEntries):
						#Calculate mean and sigma
						mean = simDataSet.get(entry).getRealValue("energyVar")
						sigma = math.sqrt(pow(testAlpha*mean,2)+pow(testBeta,2)*mean+pow(testGamma,2))
						#Generate a smeared values
						for i in range(0,valsToGen):
							smearedVal = testSlope*(numpy.random.normal(mean,sigma,1) - testOffset)
							#If the smeared value is in the import range, add to the smeared data set
							if (smearedVal >= importLowerEnergyBound) and (smearedVal <= importUpperEnergyBound):
								energyVar.setVal(smearedVal)
								smearedSimDataSet.add(argSet)
					
					#Make sim pdf
					if useSimKeysPdf==0:
						(smearedSimDataSet.get().find("energyVar")).setBins(nBins)
						smearedSimDataHist = smearedSimDataSet.binnedClone()
						simPdf = ROOT.RooHistPdf("simPdf","simPdf",argSet,smearedSimDataHist,1) #1 specifies interpolation order
					else:
						simPdf = ROOT.RooKeysPdf("simPdf","simPdf",energyVar,smearedSimDataSet,ROOT.RooKeysPdf.NoMirror,0.9)
					
					
					
					##############
					##Make Model##
					##############
					model = ROOT.RooAddPdf("model","model",ROOT.RooArgList(bgndDataPdf,simPdf),ROOT.RooArgList(scaledBgndEntriesVar,sourceCountsVar))
					model.fixCoefRange("fitRange")
					
					
					
					###########
					#Fit model#
					###########
					#Fit to sourceDataHist if using binned fits, may speed things up.
					res = model.fitTo(sourceDataSet,
										ROOT.RooFit.Extended(1),
										ROOT.RooFit.Range("fitRange"),
										ROOT.RooFit.SumCoefRange("fitRange"),
										ROOT.RooFit.NumCPU(numCPUForFits),
										ROOT.RooFit.Verbose(0),
										ROOT.RooFit.PrintLevel(0),
										ROOT.RooFit.Save(True)
									)



					########
					##Plot##
					########
					c1.cd()
					c1.SetTitle("Channel "+str(dataTreeChannelNum)+" fit "+str(stepNum)+" of "+str(totalSteps))
					#Make frame, set title
					frame = energyVar.frame(fitRangeMin,fitRangeMax,nBins)
					frame.SetTitle("alpha="+str(testAlpha)+", beta="+str(testBeta)+", gamma="+str(testGamma)+", slope="+str(testSlope)+", offset="+str(testOffset))
					#Plot source data
					sourceDataSet.plotOn(frame,ROOT.RooFit.Name("Source"),ROOT.RooFit.MarkerColor(1),ROOT.RooFit.FillColor(0) )
					#Plot model and components
					model.plotOn(frame,ROOT.RooFit.Name("Fit"),ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.FillColor(0),ROOT.RooFit.ProjWData(sourceDataSet))
					model.plotOn(frame,ROOT.RooFit.Name("Bgnd"),ROOT.RooFit.Components("bgndDataPdf"),ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.FillColor(0),ROOT.RooFit.ProjWData(sourceDataSet))
					model.plotOn(frame,ROOT.RooFit.Name("Sim"),ROOT.RooFit.Components("simPdf"),ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue),ROOT.RooFit.FillColor(0),ROOT.RooFit.ProjWData(sourceDataSet))
					#Draw
					frame.Draw()

					#Make legend
					leg = ROOT.TLegend(0.65,0.65,0.95,0.92);
					sourceDataObj = frame.findObject("Source");
					bgndDataObj = frame.findObject("Bgnd");
					simDataObj = frame.findObject("Sim");
					modelDataObj = frame.findObject("Fit")
					leg.AddEntry(sourceDataObj, "Source data","P");
					leg.AddEntry(bgndDataObj, "Bgnd data","L");
					leg.AddEntry(simDataObj, "Sim","L");
					leg.AddEntry(modelDataObj, "Bgnd+Sim","L");
					leg.Draw("same")

					#Update plots
					c1.Modified()
					c1.Update()
					ROOT.gSystem.ProcessEvents()


					#############
					#Fill TTree##
					#############			
					alpha[0]=testAlpha
					beta[0]=testBeta
					gamma[0]=testGamma
					
					slope[0]=testSlope

					offset[0]=testOffset
					
					fitNLL[0] = res.minNll();
					
					fitStatus[0] = res.status()
					print("Fit nll:"+str(fitNLL[0])+", fit status: "+str(fitStatus[0]))
					
					#Extended
					fitBgndCounts[0] = scaledBgndEntriesVar.getVal()
					fitSourceCounts[0] = sourceCountsVar.getVal()
					fitSourceCountsError[0] = sourceCountsVar.getError()
								   
					fitTree.Fill()
					
					if (fitNLL[0] < bestNLL) and fitStatus[0]==0:
						bestNLL=fitNLL[0]
						bestSlope=slope[0]
						bestOffset=offset[0]
						bestAlpha=alpha[0]
						bestBeta=beta[0]
						bestGamma=gamma[0]
						bestBgndCounts=fitBgndCounts[0]
						bestSourceCounts=fitSourceCounts[0]
						bestSourceCountsError=fitSourceCountsError[0]
						bestFitStatus=fitStatus[0]
						
						print("Found new lower NLL value: "+str(bestNLL)+", alpha="+str(bestAlpha)+", beta="+str(bestBeta)+", gamma="+str(bestGamma)+", slope="+str(bestSlope)+", offset="+str(bestOffset)+"\n\n\n")
						
						#Updated code to only write new lower NLL plots
						outputFile.cd()
						c1.Write("alpha_"+str(testAlpha)+"beta_"+str(testBeta)+"gamma_"+str(testGamma)+"_slope"+str(testSlope)+"_offset"+str(testOffset))
					
					#Memory clean-up. Probably no longer necessary since leaks were fixed.
					frame.Delete()
					del frame
					simPdf.Delete()
					del simPdf
					model.Delete()
					del model
					gc.collect()
					
print("Found lowest NLL value of "+str(bestNLL)+" for alpha="+str(bestAlpha)+", beta="+str(bestBeta)+", gamma="+str(bestGamma)+", slope="+str(bestSlope)+", offset="+str(bestOffset)+"\n\n\n")			



#####################
##Gen profile plots##
#####################
nEntries = fitTree.GetEntries()
#Expects NLL variable at end!
fitVarNames=["alpha","beta","gamma","slope","offset","nll"]
fitVars=[]

for entryNum in range(0,nEntries):
	fitTree.GetEntry(entryNum)
	fitVars.append([fitTree.alpha,fitTree.beta,fitTree.gamma,fitTree.slope,fitTree.offset,fitTree.fitNLL])

graphs=[]
#Step through every component but NLL
for parNum in range(0, len(fitVarNames)-1):

	graphs.append(ROOT.TGraph())
	graphs[-1].SetTitle(fitVarNames[parNum])
	graphs[-1].SetName(fitVarNames[parNum])
	graphs[-1].GetXaxis().SetTitle(fitVarNames[parNum])
	graphs[-1].GetYaxis().SetTitle("NLL")
	nPoints=0
	
	#Make list of unique values of this par
	uniqueElementList=[]
	for fitValue in range(0,len(fitVars)):
		if fitVars[fitValue][parNum] not in uniqueElementList:
			uniqueElementList.append(fitVars[fitValue][parNum])

	#Now calculate min nll for each value 
	minNlls=[]
	for uniqueElement in uniqueElementList:
		nlls=[]
		for fit in fitVars:
			if fit[parNum]==uniqueElement:
				nlls.append(fit[len(fitVars[0])-1])

		minNll = min(nlls)
		minNlls.append(minNll)

	#Offset so lowerest NLL occurs at zero
	nllOffset=min(minNlls)
	for uniqueElement in uniqueElementList:
		graphs[-1].SetPoint(nPoints,uniqueElement,minNlls[nPoints]-nllOffset)
		nPoints+=1

	graphs[-1].SetMarkerStyle(20)



##################
#Fill bestFitTree#
##################
alpha[0]=bestAlpha
beta[0]=bestBeta
gamma[0]=bestGamma
slope[0]=bestSlope
offset[0]=bestOffset
fitNLL[0] = bestNLL
fitBgndCounts[0] = bestBgndCounts
fitSourceCounts[0] = bestSourceCounts
fitSourceCountsError[0] = bestSourceCountsError
fitStatus[0] = bestFitStatus
bestFitTree.Fill()



#####################
#Write trees to file#
#####################
outputFile.cd()
fitTree.Write()
bestFitTree.Write()
for graph in graphs:
	graph.Write()
outputFile.Close() 
	
		
		
		
		

