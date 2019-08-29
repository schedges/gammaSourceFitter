#This version of the fitting code loops over values of different fitting parameters in discrete steps,
#fits the ampltitude of the simulation PDF, and returns the NLL value. This utilizes a binned simulation
#and bgnd PDFs. Care should be taken to test the effect of binning on parameter values.
#
#This code outputs three things to the specified output file
#	fitTree: containing fit parameters for each fit
#	bestFitTree: containing fit parameters for fit with lowest NLL only
#	canvases of fits: Currently only saves fits it finds with a lower NLL than the current lowest NLL value,
#		so last canvas saved should correspond to the best fit.
import ROOT
import array
import numpy
import math
import gc

#Specify number of bins for PDFs, over the entire import range (lowerEnergyBound, upperEnergyBound)
nBins=200

#Specify how many threads to use for fit
numCPUForFits=1

##################################
##Specify data format/names here##
##################################
#For data trees with multiple channels, specify the channel number here
dataTreeChannelNum=0

#Specify name of bgnd file and source file (and path if not in current directory).
#This code assumes the bgnd data set and source data set have the same tree name, branches, etc.
bgndFilename = "SIS3316Raw_20180823033124.root"
sourceFilename = "SIS3316Raw_20180821164027.root"

#Specify name of data tree
dataTreeName = "BDTree"
#Specify name and type of branch corresponding to energy
dataTreeEnergyBranchName ="integral"
dataTreeEnergyBranchType = 'd'
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
simFilename = "dumn1-naiqf-662_v1-1e8.root"

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
#Specify the range to fit over (in terms of dataTreeEnergy).
#NOTE: you probably want the minimum value to be above the threshold region.
fitRangeMin=2000
fitRangeMax=14000
#Specify range for importing energy variables from both dataTrees and simTrees.
#NOTE: this should be larger than the fit range, as data from outside the fit range
#can be brought into it when applying smearing. I've also had problems when the bounds
#includes 0, so I start just above that.
#NOTE: These values are inclusive
lowerEnergyBound = 1
upperEnergyBound = 20000


#############################
##keV->ADC Start Parameters##
##c3###########################
#Typical parameterization of light output as a fn of energy, L(E), above 40 keV, is
#L = c * (E-E_0)
#where c=slope, the conversion factor, and E_0=offset,  a small offset,
#usually between 0 and 25 keV

#Specify slope range and number of steps to include between that range here
slopeMin=15
slopeMax=21
numSlopeSteps=4

#Specify offset here. Note that because of the way I define offset later, I think this
#should be a negative keV value
offsetMin=-25
offsetMax=0
numOffsetSteps=2


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
alphaMin = 0.1
alphaMax = 0.15
numAlphaSteps = 3

betaMin = 0.05
betaMax = 0.15
numBetaSteps = 5

gammaMin = 0.0
gammaMax = 0.065
numGammaSteps= 3

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


###################
##RooFit Settings##
###################
#Reduces warning messages. I suggest leaving these on at first to make sure there's no serious
#warning/error messages
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

#Controls the number of steps used for integration. I found fits failing with the default
#so I increased it to 30
ROOT.RooAbsReal.defaultIntegratorConfig().getConfigSection("RooIntegrator1D").setRealValue("maxSteps",30);

#Set integration step sizes, useful for balancing speed vs. accuracy
ROOT.RooAbsReal.defaultIntegratorConfig().setEpsAbs(1e-6)
ROOT.RooAbsReal.defaultIntegratorConfig().setEpsRel(1e-6)


################
##For plotting##
################
c1 = ROOT.TCanvas("c1","Channel "+str(dataTreeChannelNum)+" fits")


#######################MAIN CODE STARTS HERE##########################


################################
##Set up observable, fit range##
################################
energyVar = ROOT.RooRealVar("energyVar","energyVar",lowerEnergyBound,upperEnergyBound)
energyVar.setRange("fitRange",fitRangeMin,fitRangeMax)
binning = ROOT.RooBinning(nBins,lowerEnergyBound,upperEnergyBound,"binning")
energyVar.setBinning(binning)
#We'll apply this cut to the source and bgnd data sets solely for determining the expected
#counts in those PDFs for use in the fit
cut="energyVar>="+str(fitRangeMin)+"&&energyVar<="+str(fitRangeMax)


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
nBgndEntries=bgndTree.GetEntries()
print("Found "+str(nBgndEntries)+" in bgnd tree")

#Step through bgnd tree, adding data if meets energy and channel criteria
for entry in range(0,nBgndEntries):
	bgndTree.GetEntry(entry)
	#Check if this belongs to the channel number we're fitting data for
	if dataTreeChannel[0] == dataTreeChannelNum:
		#Check if the value is in the import range
		if (dataTreeEnergy[0] >= lowerEnergyBound) and (dataTreeEnergy[0] <= upperEnergyBound):
			#Set value of observable
			energyVar.setVal(dataTreeEnergy[0])
			#Add to data set
			bgndDataSet.add(ROOT.RooArgSet(energyVar))

#Get length or run for normalization
#Get first entry
bgndTree.GetEntry(0)
startTime = dataTreeTime[0]
#Get last entry
bgndTree.GetEntry(nBgndEntries-1)
endTime = dataTreeTime[0]
#Calculate run length
bgndLengthOfRun = endTime-startTime
print("Bgnd run length: "+str(bgndLengthOfRun))

#Get number of counts in fit range
reducedDataSet = bgndDataSet.reduce(ROOT.RooFit.Cut(cut))
bgndCountsInFitRange = reducedDataSet.numEntries()
print("Found "+str(bgndCountsInFitRange)+" bgnd entries in the fit range")

#Make background pdf
(bgndDataSet.get().find("energyVar")).setBins(nBins)
bgndDataHist = bgndDataSet.binnedClone()
bgndDataHistPdf = ROOT.RooHistPdf("bgndDataHistPdf","bgndDataHistPdf",ROOT.RooArgSet(energyVar),bgndDataHist,1) #1 specifies interpolation order


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
sourceDataSet = ROOT.RooDataSet("sourceDataSet","sourceDataSet",ROOT.RooArgSet(energyVar))

#Get number of entries
nSourceEntries=sourceTree.GetEntries()
print("Found "+str(nSourceEntries)+" entries in source data set")

#Step through source data set, adding entries if meet the appropriate energy and channel criteria
for entry in range(0,nSourceEntries):
	sourceTree.GetEntry(entry)
	if dataTreeChannel[0] == dataTreeChannelNum:
		if (dataTreeEnergy[0] >= lowerEnergyBound) and (dataTreeEnergy[0] <= upperEnergyBound):
			energyVar.setVal(dataTreeEnergy[0])
			sourceDataSet.add(ROOT.RooArgSet(energyVar))

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
print("Found "+str(sourceCountsInFitRange)+" source entries in the fit range")

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
simDataSet = ROOT.RooDataSet("simDataSet","simDataSet",ROOT.RooArgSet(energyVar))

#Get number of entries
nSimEntries=simTree.GetEntries()
print("Found "+str(nSimEntries)+" entries in simulation tree")

#Step through, loading entries to data set if meet energy and channel criteria
for entry in range(0,nSimEntries):
	simTree.GetEntry(entry)
	if simTreeChannel[0]==simTreeChannelNum:
		if (simTreeEnergy[0] >= lowerEnergyBound) and (simTreeEnergy[0] <= upperEnergyBound):
			energyVar.setVal(simTreeEnergy[0])
			simDataSet.add(ROOT.RooArgSet(energyVar))

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
alphaError=array.array('d',[0])
beta=array.array('d',[0])
betaError=array.array('d',[0])
gamma=array.array('d',[0])
gammaError=array.array('d',[0])
slope=array.array('d',[0])
slopeError=array.array('d',[0])
offset=array.array('d',[0])
offsetError=array.array('d',[0])
fitNLL=array.array('d',[0])
fitBgndCounts=array.array('d',[0])
fitSourceCounts=array.array('d',[0])
fitSourceCountsError=array.array('d',[0])
fitStatus=array.array('i',[0])

fitTree.Branch('alpha',alpha,'alpha/D')
fitTree.Branch('alphaError',alphaError,'alphaError/D')
fitTree.Branch('beta',beta,'beta/D')
fitTree.Branch('betaError',betaError,'betaError/D')
fitTree.Branch('gamma',gamma,'gamma/D')
fitTree.Branch('gammaError',gammaError,'gammaError/D')
fitTree.Branch('slope',slope,'slope/D')
fitTree.Branch('slopeError',slopeError,'slopeError/D')
fitTree.Branch('offset',offset,'offset/D')
fitTree.Branch('offsetError',offsetError,'offsetError/D')
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
bestFitTree.Branch('alphaError',alphaError,'alphaError/D')
bestFitTree.Branch('beta',beta,'beta/D')
bestFitTree.Branch('betaError',betaError,'betaError/D')
bestFitTree.Branch('gamma',gamma,'gamma/D')
bestFitTree.Branch('gammaError',gammaError,'gammaError/D')
bestFitTree.Branch('slope',slope,'slope/D')
bestFitTree.Branch('slopeError',slopeError,'slopeError/D')
bestFitTree.Branch('offset',offset,'offset/D')
bestFitTree.Branch('offsetError',offsetError,'offsetError/D')
bestFitTree.Branch('fitNLL',fitNLL,'fitNLL/D')
bestFitTree.Branch('fitBgndCounts',fitBgndCounts,'fitBgndCounts/D')
bestFitTree.Branch('fitSourceCounts',fitSourceCounts,'fitSourceCounts/D')
bestFitTree.Branch('fitSourceCountsError',fitSourceCountsError,'fitSourceCountsError/D')
bestFitTree.Branch('fitStatus',fitStatus,'fitStatus/I')

bestNLL=1000000 #Set to something large so hopefully first fit has a lower nll and we plot it
bestSlope=0
bestSlopeError=0
bestOffset=0
bestOffsetError=0
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
					#Assume we're within 20% of the true values
					sourceCountsVar = ROOT.RooRealVar("sourceCountsVar","sourceCountsVar",expectedSourceCounts,0.8*expectedSourceCounts,1.2*expectedSourceCounts)
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
						sigma = 1./2.355 * math.sqrt(pow(testAlpha*mean,2)+pow(testBeta,2)*mean+pow(testGamma,2))
						#Generate a smeared values
						smearedVal = testSlope*(numpy.random.normal(mean,sigma,1) - testOffset)
						#If the smeared value is in the import range, add to the smeared data set
						if (smearedVal >= lowerEnergyBound) and (smearedVal <= upperEnergyBound):
							energyVar.setVal(smearedVal)
							smearedSimDataSet.add(ROOT.RooArgSet(energyVar))
					

					#Make background pdf
					(smearedSimDataSet.get().find("energyVar")).setBins(nBins)
					smearedSimDataHist = smearedSimDataSet.binnedClone()
					simPdf = ROOT.RooHistPdf("simPdf","simPdf",ROOT.RooArgSet(energyVar),smearedSimDataHist,1) #1 specifies interpolation order

					
					
					##############
					##Make Model##
					##############
					model = ROOT.RooAddPdf("model","model",ROOT.RooArgList(bgndDataHistPdf,simPdf),ROOT.RooArgList(scaledBgndEntriesVar,sourceCountsVar))
					model.fixCoefRange("fitRange")
					
					###########
					#Fit model#
					###########
					res = model.fitTo(sourceDataSet,
									  ROOT.RooFit.Extended(1),
									  ROOT.RooFit.Range("fitRange"),
									  ROOT.RooFit.SumCoefRange("fitRange"),
									  ROOT.RooFit.NumCPU(numCPUForFits),
									  ROOT.RooFit.Verbose(0),
									  ROOT.RooFit.PrintLevel(-1),
									  ROOT.RooFit.Save(True)
									)

					########
					##Plot##
					########
					c1.cd()
					#Make frame, set title
					frame = energyVar.frame(fitRangeMin,fitRangeMax,nBins)
					frame.SetTitle("alpha="+str(testAlpha)+", beta="+str(testBeta)+", gamma="+str(testGamma)+", slope="+str(testSlope)+", offset="+str(testOffset))
					#Plot source data
					sourceDataSet.plotOn(frame,ROOT.RooFit.Name("Source"),ROOT.RooFit.MarkerColor(1),ROOT.RooFit.FillColor(0) )
					#Plot model and components
					model.plotOn(frame,ROOT.RooFit.Name("Fit"),ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.FillColor(0),ROOT.RooFit.ProjWData(sourceDataSet))
					model.plotOn(frame,ROOT.RooFit.Name("Bgnd"),ROOT.RooFit.Components("bgndDataHistPdf"),ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.FillColor(0),ROOT.RooFit.ProjWData(sourceDataSet))
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
					alphaError[0]=(alphaMax-alphaMin)/numAlphaSteps
					beta[0]=testBeta
					betaError[0]=(betaMax-betaMin)/numBetaSteps
					gamma[0]=testGamma
					gammaError[0]=(gammaMax-gammaMin)/numBetaSteps
					
					slope[0]=testSlope
					slopeError[0]=(slopeMax-slopeMin)/numSlopeSteps

					offset[0]=testOffset
					offsetError[0]=(offsetMax-offsetMin)/numOffsetSteps
					
					fitNLL[0] = res.minNll();
					print("Fit nll:"+str(res.minNll()))
					
					#Extended
					fitBgndCounts[0] = scaledBgndEntriesVar.getVal()
					fitSourceCounts[0] = sourceCountsVar.getVal()
					fitSourceCountsError[0] = sourceCountsVar.getError()
							     
					fitStatus[0] = res.status()
							      
					fitTree.Fill()
					
					if (fitNLL[0] < bestNLL):
						bestNLL=fitNLL[0]
						bestSlope=slope[0]
						bestSlopeError=slopeError[0]
						bestOffset=offset[0]
						bestOffsetError=offsetError[0]
						bestAlpha=alpha[0]
						bestAlphaError=alphaError[0]
						bestBeta=beta[0]
						bestBetaError=betaError[0]
						bestGamma=gamma[0]
						bestGammaError=gammaError[0]
						bestBgndCounts=fitBgndCounts[0]
						bestSourceCounts=fitSourceCounts[0]
						bestSourceCountsError=fitSourceCountsError[0]
						bestFitStatus=fitStatus[0]
						
						print("Found new lower NLL value: "+str(bestNLL)+", alpha="+str(bestAlpha)+", beta="+str(bestBeta)+", gamma="+str(bestGamma)+", slope="+str(bestSlope)+", offset="+str(bestOffset)+"\n\n\n")
						
						#Updated code to only write new lower NLL plots
						outputFile.cd()
						c1.Write("alpha_"+str(testAlpha)+"beta_"+str(testBeta)+"gamma_"+str(testGamma)+"_slope"+str(testSlope)+"_offset"+str(testOffset))
					
					
					#There's a memory leak, trying to fix
					frame.Delete()
					del frame
					simPdf.Delete()
					del simPdf
					model.Delete()
					del model
					gc.collect()
					
print("Found lowest NLL value of "+str(bestNLL)+" for alpha="+str(bestAlpha)+", beta="+str(bestBeta)+", gamma="+str(bestGamma)+", slope="+str(bestSlope)+", offset="+str(bestOffset)+"\n\n\n")			

##################
#Fill bestFitTree#
##################
alpha[0]=bestAlpha
alphaError[0]=bestAlphaError
beta[0]=bestBeta
betaError[0]=bestBetaError
gamma[0]=bestGamma
gammaError[0]=bestGammaError
slope[0]=bestSlope
slopeError[0]=bestSlopeError
offset[0]=bestOffset
offsetError[0]=bestOffsetError 
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
outputFile.Close() 
		
		
		
		
		

