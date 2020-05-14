#!/var/phy/project/phil/sch52/nubeCalibrations/env/bin/python3.6
# I'm running using the above virtual environment, and by specifying it and 
# making your python script executable, you can run without needing to activate
# it.
#
# Description: Emcee and RooFit based code to calibrated liquid scintillators
# using source(s), a background spectrum, and simulation file(s).
#
#Required packages that need to be installed (via pip):
#    emcee
#    corner 
#    matplotlib 
#    tqdm
#
#Usage:
# ./fitSim_emcee_simultaneous <path to folder> <channel num>
#
# Output plots and sampler will be written to <path to folder>. Assumes 
# <path to folder> contains a sub folder named ls<channel num>, and inside 
# ls<channel num> should be subfolders with each source name and a root file 
# named "bgnd.root". Inside each source subfolder should be a file named
# <source name>.root and <source name>sim.root. 
#
# Notes: 
#   - sources to calibrate in ls<channel num> folder are hardcoded, as is the 
#     range to fit.
#   - I've had slowdowns with RooFit when calling createNll too many times. This
#     appears to be a known bug. The solution is to take advantage of the emcee
#     h5py backend and break up the run into multiple
#   - If your normalization is off, you may want to let the background pdf amp
#     float.

import ROOT
import numpy
import array
import matplotlib
import math
import gc
import csv as csvlib
import sys
import emcee
import corner
import os
from matplotlib import pyplot as plt
from multiprocessing import Pool

#Make sure process is single-threaded, we'll use multiprocessing library instead
os.environ["OMP_NUM_THREADS"] = "1"

##################
##Get user input##
##################
calibrationFolder=sys.argv[1]
if not calibrationFolder.endswith('/'):
  calibrationFolder+="/"
#Specify where to write plots
plotPath=calibrationFolder

calibrationChannel = sys.argv[2]

calibrationSources=["na22","cs137","co60"]

# Specify data tree name (Mine changes depending on what channel we're fitting), 
# and branch names/types containing channel, energy quantitiy, and run time
leadTreeChannels=[1,2,12,13]
ironTreeChannels=[0,3,14,15]
if int(calibrationChannel) in leadTreeChannels:
  dataTreeName="leadTree"
else:
  dataTreeName="ironTree"
dataTreeChannelBranchName = "channelID"
dataTreeChannelBranchType='i'
dataTreeEnergyBranchName ="integral"
dataTreeEnergyBranchType = 'd'
dataTreeTimeBranchName = "waveformStartTimeSecondsPart"
dataTreeTimeBranchType = 'd'

#Specify sim tree name, and branch name/type of energy quantity
simTreeName = "combinedEnergyTree"
simTreeEnergyBranchName = "energy"
simTreeEnergyBranchType = 'd'

################
##emcee set-up##
################
ndim=5 #Alpha, Beta, Gamma, Slope, Offset
#Add one dimension for each source to represent amplitude
for source in calibrationSources:
  ndim+=1
  
#Parameter names for plots
labels=['alpha','beta','gamma','slope','offset'] 
for source in calibrationSources:
  labels.append("amp_"+source)

#Minimum values for the above parameters
mins=[0.0,0.05,0.0,45,-5.]
maxes=[0.09,9.0,5.0,85,5.]
#I'm limiting amplitudes to 0.9-1.01 expected amplitude
for source in calibrationSources:
  mins.append(0.90)
  maxes.append(1.01)
  
nwalkers=500 #Should be at least 2x ndim, but emcee recommends "hundreds"
#Ideally several auto-correlation length
nBurnInSteps=100
#Ideally 10x autocorrelation length
nSteps=1000

#######################
#Data loading settings#
####################### 
# Set to one to limit the number of source data entries to 
# maxSourceDataEntriesToUse. I think ideally there should be a similar number of
# entries in each source run so the source fits are weighted ~equally. Also 
# decreasing this speeds up fits.
limitToMaxEntries=1
maxSourceDataEntriesToUse=160000

#Same thing for sims, this really just affects the speed 
limitToMaxSimEntries=1
maxSimEntriesToUse=100000

# How many smeared samples to generate for each real simulation sample. Smooths 
# PDF but increases fit time 
valsToGen=20

# Either "keys" or "binned". 
# - "keys" makes a RooKeysPdf of the background spectrum, which ensures there 
# are no zero bins but I've noticed has problems near rapidly changing regions 
# such as near threshold. 
# - "binned" models the background as a RooHistPdf, and will have issues if 
# stats are poor or binning is too fine and there are bins with zero-content.
bgndPdfType="keys"

#Specify binning for source and sim pdfs, and bgnd pdf if above is "binned"
nBins=200

####################################
##Specify data import / fit ranges##
####################################
# Specify the range to fit over (in terms of data energy quantity).
# NOTE: you probably want the minimum value to be above the threshold region.
fitRangeMin=12000
fitRangeMax=150000
# Specify range for importing energy variables from both dataTrees and simTrees.
# NOTE: this should be larger than the fit range, as data from outside the fit 
# range can be brought into it when applying smearing. I've also had problems 
# when the bounds include 0, so I start just above that.
# NOTE: These values are inclusive
lowerEnergyBound = 1
upperEnergyBound = 200000



################################################################################
##Assuming you followed the specified folder structure, you shouldn't need to###
##modify anything below here...#################################################
################################################################################
#* other than changing keyspdf setting or interpolation options



###################
##RooFit Settings##
###################
# Reduces warning messages.
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

# Controls the number of steps used for integration. I needed to increase this
# in some fits when generating the NLL
ROOT.RooAbsReal.defaultIntegratorConfig().getConfigSection("RooIntegrator1D").setRealValue("maxSteps",30);

# Set integration step sizes, useful for balancing speed vs. accuracy
ROOT.RooAbsReal.defaultIntegratorConfig().setEpsAbs(1e-7)
ROOT.RooAbsReal.defaultIntegratorConfig().setEpsRel(1e-7)

##################################
##Specify data format/names here##
##################################
dataTreeChannelNum=int(calibrationChannel)

basePath=calibrationFolder+"ls"+calibrationChannel
bgndFilename = basePath+"/bgnd.root"
print("Bgnd file is: "+bgndFilename)
baseSourcePaths=[]
sourceFilenames=[]
simFilenames=[]
for source in calibrationSources:
  baseSourcePaths.append(basePath+"/"+source)
  sourceFilenames.append(baseSourcePaths[-1]+"/"+source+".root")
  print("Source file is: "+sourceFilenames[-1])
  simFilenames.append(baseSourcePaths[-1]+"/"+source+"sim.root")
  print("Sim file is: "+sourceFilenames[-1])

dataTreeEnergy = array.array(dataTreeEnergyBranchType,[0])
dataTreeChannel = array.array(dataTreeChannelBranchType,[0])
dataTreeTime = array.array(dataTreeTimeBranchType,[0])

################################
##Set up observable, fit range##
################################
energyVar = ROOT.RooRealVar("energyVar","energyVar",lowerEnergyBound,upperEnergyBound)
binning = ROOT.RooBinning(nBins,lowerEnergyBound,upperEnergyBound,"binning")
energyVar.setBinning(binning)
energyVar.setRange("fitRange",fitRangeMin,fitRangeMax)
# Make RooArgSet, can leave to memory leaks if made in a loop
argSet=ROOT.RooArgSet(energyVar)
argList=ROOT.RooArgList(energyVar)
# We'll apply this cut to source and bgnd data sets to determine the expected
# counts in those PDFs for use in the fit
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
bgndDataSet = ROOT.RooDataSet("bgndDataSet","bgndDataSet",argSet)

#Get number of entries
nBgndEntries=bgndTree.GetEntries()
print("\nFound "+str(nBgndEntries)+" in bgnd tree")

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
reducedBgndDataSet=bgndDataSet.reduce(ROOT.RooFit.Cut(cut))
bgndCountsInFitRange = reducedBgndDataSet.numEntries()
print("Found "+str(bgndCountsInFitRange)+" bgnd entries in the fit range")

#Make background pdf
if bgndPdfType=="binned":
  print("Making binned background pdf\n")
  (bgndDataSet.get().find("energyVar")).setBins(nBins)
  bgndDataHist = bgndDataSet.binnedClone()
  bgndDataPdf = ROOT.RooHistPdf("bgndDataPdf","bgndDataPdf",argSet,bgndDataHist,1) #1 specifies interpolation order
else:
  print("Making RooKeysPdf background pdf\n")
  bgndDataPdf = ROOT.RooKeysPdf("bgndDataPdf","bgndDataPdf",energyVar,bgndDataSet,ROOT.RooKeysPdf.NoMirror,1.2)

########################
##Load source data set## 
########################
#Open source file
sourceDataSets=[]
sourceDataHists=[]
sourceCountsInFitRanges=[]
sourceLengthOfRuns=[]
srcDataHistPdfs=[]
for sourceNum in range(0,len(calibrationSources)):

  #Get name of source 
  source=calibrationSources[sourceNum]
  
  #Open file 
  sourceFile = ROOT.TFile(sourceFilenames[sourceNum],"READ")
  
  #Get source tree
  sourceTree = sourceFile.Get(dataTreeName)
  
  #Set branches
  sourceTree.SetBranchAddress(dataTreeEnergyBranchName,dataTreeEnergy)
  sourceTree.SetBranchAddress(dataTreeTimeBranchName,dataTreeTime)
  sourceTree.SetBranchAddress(dataTreeChannelBranchName,dataTreeChannel)

  #Create RooDataSet for holding source data
  sourceDataSets.append(ROOT.RooDataSet("sourceDataSet_"+source,"sourceDataSet_"+source,argSet))

  #Get number of entries
  nSourceEntries=sourceTree.GetEntries()
  print("Found "+str(nSourceEntries)+" entries in "+source+" data set")

  #Step through source data set, adding entries if meet the appropriate energy and channel criteria
  if limitToMaxEntries==1:
    if nSourceEntries>maxSourceDataEntriesToUse:
      nSourceEntries=maxSourceDataEntriesToUse
      print("Limiting max "+source+" data entries to "+str(nSourceEntries))
  for entry in range(0,nSourceEntries):
    sourceTree.GetEntry(entry)
    if dataTreeChannel[0] == dataTreeChannelNum:
      if (dataTreeEnergy[0] >= lowerEnergyBound) and (dataTreeEnergy[0] <= upperEnergyBound):
        energyVar.setVal(dataTreeEnergy[0])
        sourceDataSets[sourceNum].add(argSet)

  #Get length or run for normalization
  #Get first entry
  sourceTree.GetEntry(0)
  startTime = dataTreeTime[0]
  #Get last entry
  sourceTree.GetEntry(nSourceEntries-1)
  endTime = dataTreeTime[0]
  #Calculate length of run
  sourceLengthOfRuns.append(endTime-startTime)

  #Get number of counts in fit range
  sourceCountsInFitRanges.append(sourceDataSets[sourceNum].numEntries())
  print("Found "+str(sourceCountsInFitRanges[-1])+" entries for "+source+" source the fit range\n")

  (sourceDataSets[sourceNum].get().find("energyVar")).setBins(nBins)
  sourceDataHists.append(sourceDataSets[sourceNum].binnedClone())

#################################################
##Normalize bgnd counts to length of source run##
#################################################
print("Length of bgnd run: "+str(bgndLengthOfRun)+"\n")

expectedSourceCounts=[]
scaledBgndEntriesVars=[]
for sourceNum in range(0,len(calibrationSources)):
  source=calibrationSources[sourceNum]
  
  #Scale length of background run to match source run
  scaleFactor = float(sourceLengthOfRuns[sourceNum])/float(bgndLengthOfRun)
  scaledBgndEntries = bgndCountsInFitRange*scaleFactor
  #Make a RooRealVar to hold this, and set expected bgnd counts constant
  scaledBgndEntriesVars.append(ROOT.RooRealVar("scaledBgndEntriesVar_"+source,"scaledBgndEntriesVar_"+source,scaledBgndEntries))
  scaledBgndEntriesVars[-1].setConstant()
  print("Length of "+source+" run: "+str(sourceLengthOfRuns[sourceNum]))
  print(source+" Scale factor: "+str(scaleFactor))
  print(source+" Scaled bgnd counts in fit range: "+str(scaledBgndEntries))
  print("Expected "+source+" counts in fit range: "+str(sourceCountsInFitRanges[sourceNum]-scaledBgndEntries)+"\n")
  expectedSourceCounts.append(sourceCountsInFitRanges[sourceNum]-scaledBgndEntries)

#####################
##Load Sim Data Set##
#####################
simTreeEnergy = array.array(simTreeEnergyBranchType,[0])
simDataSets=[]
simArrays=[]
for sourceNum in range(0,len(calibrationSources)):

  source=calibrationSources[sourceNum]

  #Open sim file
  simFile = ROOT.TFile(simFilenames[sourceNum],"READ")
  #Get simTree from file
  simTree = simFile.Get(simTreeName)
  #Load branches
  simTree.SetBranchAddress(simTreeEnergyBranchName,simTreeEnergy)

  #Create data set for sim
  simDataSets.append(ROOT.RooDataSet("simDataSet_"+source,"simDataSet_"+source,argSet))

  #Get number of entries
  nSimEntries=simTree.GetEntries()
  print("Found "+str(nSimEntries)+" entries in "+source+" simulation tree")
  
  if limitToMaxSimEntries==1:
    if nSimEntries>maxSimEntriesToUse:
      nSimEntries=maxSimEntriesToUse
      print("Limiting "+source+" sim to "+str(nSimEntries)+" entries\n")

  #Numpy way
  simList=[]
  for entry in range(0,nSimEntries):
    simTree.GetEntry(entry)
    simList.append(simTreeEnergy[0])
  simArrays.append(numpy.array(simList))


################
##EMCEE Set-up##
################
PLOT=0

#Make into numpy arrays 
pos_min = numpy.array(mins)
pos_max = numpy.array(maxes)

#Make ranges for corner plots 
ranges=[]
for i in range(0,len(mins)):
  entry=[]
  entry.append(mins[i])
  entry.append(maxes[i])
  ranges.append(entry)

#size of each parameter space 
psize = pos_max - pos_min 

#Generate random values within that space for each walker
pos = [pos_min + psize*numpy.random.rand(ndim) for i in range(nwalkers)]

######################EMCEE FUNCTIONS DEFINED HERE##############################
#Efficient way to apply resolution to an array of values 
f = lambda x,alpha,beta,gamma: numpy.sqrt( numpy.power(alpha*x,2) + numpy.power(beta,2)*x + numpy.power(gamma,2) )

#Returns 0 if all parameters are in their allowed range, otherwise -infinity
def lnprior(theta):
  #Check if all parameters in range 
  allParsInRange=1
  for i in range(0,len(theta)):
    if not mins[i]<theta[i]<maxes[i]:
      allParsInRange=0
  
  #Return 0.0 if all in range, otherwise return -infinity
  if allParsInRange==1:
    return 0.0
  else:
    return -numpy.inf
  
#This is our RooFit function that returns the POSITIVE log likelihood
def lnlike(theta):
  #energyVar gets modified in loop. So far hasn't caused any issues with multiprocessing
  global energyVar

  # Reset nllVal to zero because we'll add to this for each source we generate
  # an nll value for
  nllVal=0
  
  #Load theta values into parameters
  alpha=theta[0]
  beta=theta[1]
  gamma=theta[2]
  slope=theta[3]
  offset=theta[4]
  
  # Load remaining theta values for each source, generate smeared pdf, calculate 
  # nll value
  for sourceNum in range(0,len(calibrationSources)):
    
    amplitude=theta[5+sourceNum]

    #Make ampltitude var 
    sourceCountsVar = ROOT.RooRealVar("sourceCountsVar","sourceCountsVar",amplitude*expectedSourceCounts[sourceNum])
    sourceCountsVar.setConstant(1)
    
    #Make array of sigmas from res params and sim values
    sigmas=f(simArrays[sourceNum],alpha,beta,gamma)
    
    #Generate valsToGen random values for every entry in sim
    smearedArrayList=[]
    for i in range(0,valsToGen):
      smearedArrayList.append(slope*(numpy.random.normal(simArrays[sourceNum],sigmas)-offset))
    
    #Make it a flat list 
    smearedArrayArray=numpy.array(smearedArrayList)
    flatArray=smearedArrayArray.flatten()
    
    #Make smeared data set
    smearedSimDataSet=ROOT.RooDataSet("smearedSimDataSet","smearedSimDataSet",argSet)
    
    #Numpy array->TH1->RooDataHist->RooHistPdf
    #~0.03 seconds, much faster than iterating through array to fill
    w=numpy.full(flatArray.size,1.)
    h=ROOT.TH1D("h","",nBins,lowerEnergyBound,upperEnergyBound)
    h.FillN(flatArray.size,flatArray,w)
    smearedSimDataHist=ROOT.RooDataHist("smearedSimDataHist","smearedSimDataHist",argList,h)
    simPdf = ROOT.RooHistPdf("simPdf","simPdf",argSet,smearedSimDataHist,0) #1 specifies interpolation order
    h.Delete()
    del h
    
    ##Make Model
    pdfList = ROOT.RooArgList(bgndDataPdf,simPdf)
    ampList = ROOT.RooArgList(scaledBgndEntriesVars[sourceNum],sourceCountsVar)
    model = ROOT.RooAddPdf("model","model",pdfList,ampList)
    model.fixCoefRange("fitRange")
    
    #Compute nll
    nll = model.createNLL(sourceDataHists[sourceNum],
      ROOT.RooFit.Extended(1),
      ROOT.RooFit.Verbose(0),
      ROOT.RooFit.Range("fitRange"),
      ROOT.RooFit.SumCoefRange("fitRange"),
      ROOT.RooFit.NumCPU(1)
    )

    #Make NLL positive
    nllVal += (-1*nll.getVal())
    
    if (PLOT==1):
      try:
        c1
      except NameError:
        c1=ROOT.TCanvas("c1","c1")
      
      #Reduce integrator for plotting, massively speeds plotting
      ROOT.RooAbsReal.defaultIntegratorConfig().setEpsAbs(1e-5)
      ROOT.RooAbsReal.defaultIntegratorConfig().setEpsRel(1e-5)
      
      frame = energyVar.frame(lowerEnergyBound,upperEnergyBound,nBins)
      frame.SetTitle(calibrationSources[sourceNum]+": alpha="+str(alpha)+", beta="+str(beta)+", gamma="+str(gamma)+", slope="+str(slope)+", offset="+str(offset))
      
      #Plot source data
      sourceDataSets[sourceNum].plotOn(
        frame,
        ROOT.RooFit.Name("Source"),
        ROOT.RooFit.MarkerColor(1),
        ROOT.RooFit.FillColor(0)
      )
      #Plot components
      model.plotOn(
        frame,
        ROOT.RooFit.Name("Bgnd"),
        ROOT.RooFit.Components("bgndDataPdf"),
        ROOT.RooFit.LineColor(ROOT.kSolid),
        ROOT.RooFit.FillColor(0),
        ROOT.RooFit.ProjWData(sourceDataSets[sourceNum])
      )
      model.plotOn(
        frame,ROOT.RooFit.Name("Sim"),
        ROOT.RooFit.Components("simPdf"),
        ROOT.RooFit.LineColor(ROOT.kRed),
        ROOT.RooFit.FillColor(0),
        ROOT.RooFit.ProjWData(sourceDataSets[sourceNum]),
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
      c1.SaveAs(plotPath+"bestFit_simultaneous_"+calibrationSources[sourceNum]+"_ch"+calibrationChannel+".pdf")
      
      #Reset integrator for step size
      ROOT.RooAbsReal.defaultIntegratorConfig().setEpsAbs(1e-7)
      ROOT.RooAbsReal.defaultIntegratorConfig().setEpsRel(1e-7)
      
      #Memory management for plotting
      frame.Delete()
      del Frame()
      
    #Memory management
    smearedSimDataSet.reset()
    smearedSimDataSet.Delete()
    del smearedSimDataSet
    pdfList.Delete()
    del pdfList
    ampList.Delete()
    del ampList
    sourceCountsVar.Delete()
    del sourceCountsVar
    smearedSimDataHist.Delete()
    del smearedSimDataHist
    simPdf.Delete()
    del simPdf
    model.Delete()
    del model
    nll.Delete()
    del nll
    gc.collect()
    
    
  #Return total nllval from all sources
  return nllVal
  
#Calls out roofit function, makes sure output is not infinite and parameters in allowed range 
def lnprob(theta):
  lp = lnprior(theta)
  if not numpy.isfinite(lp):
    return -numpy.inf
  return lp + lnlike(theta)



##########################MC RUN STARTS HERE####################################

#Single-threaded
#sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)
#print("Starting burn in...")
#pos, prob, state  = sampler.run_mcmc(pos, nBurnInSteps, progress=True)
#sampler.reset()
#print("Burn-in complete!")
#pos, prob, state  = sampler.run_mcmc(pos, nSteps, progress=True)

#Parallel processing
with Pool() as pool:
  sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,pool=pool)
  #Burn in
  print("Starting burn in...")
  pos, prob, state  = sampler.run_mcmc(pos, nBurnInSteps, progress=True)
  print("Burn-in complete! Mean acceptance fraction: {0:.3f}".format(numpy.mean(sampler.acceptance_fraction)))
  sampler.reset()
  pos, prob, state  = sampler.run_mcmc(pos,nSteps,progress=True)

###########################MC PLOTTING HERE#####################################
#My computer doesn't have python-tk, so I can't view plots and had to save them
#as a PDF to look at the results 
matplotlib.use('PDF')

#GET THE LL VALUES--From grayson's code--do this first in case plotting fails.
samples=sampler.flatchain
lnprobs = sampler.lnprobability[:,:]
flatLnprobs = lnprobs.reshape(-1)
with open(plotPath+"sampler_simultaneous_ch"+calibrationChannel+".csv", 'w') as sampleOutFile:
  theWriter = csvlib.writer(sampleOutFile, delimiter=',')
  for sampleLine, llval in zip(samples, flatLnprobs):
    theWriter.writerow(numpy.append(sampleLine,llval))

#MAKE TRACE PLOTS OF EACH PARAMATER
fig = plt.figure(figsize=(10,ndim*2))
gs = fig.add_gridspec(ndim,1)
plt.subplots_adjust(hspace=0.4)
for i in range(0,ndim):
  axes = fig.add_subplot(gs[i,:])
  axes.plot(sampler.chain[:,:,i].T, '-', color='k', alpha=0.3)
  axes.set_title(labels[i])
plt.savefig(plotPath+"traceplots_simultaneous_ch"+calibrationChannel+".pdf")

#CORNER PLOT HERE 
samples=sampler.flatchain
fig = corner.corner(samples, labels=labels, ranges=ranges, quantiles=[0.16,0.5,0.84],show_titles=True,title_kwargs={'fontsize':12})
fig.savefig(plotPath+"corner_simultaneous_ch"+calibrationChannel+".pdf")

#CALCULATE QUANTILES HERE
bestFitValues=[]
for i in range(ndim):
  mcmc = numpy.percentile(sampler.chain[:,:,i],[16, 50, 84])
  q = numpy.diff(mcmc)
  print(labels[i]+": "+str(mcmc[1])+"+"+str(q[0])+" -"+str(q[1])+"\n")
  bestFitValues.append(mcmc[1])

#Plot best fit 
PLOT=1
lnlike(bestFitValues) #Unpack list to arguments for our call to lnlike

#Print out stats
print("Mean acceptance fraction: {0:.3f}".format(numpy.mean(sampler.acceptance_fraction))+"\n")
print("Mean autocorrelation time: {0:.3f} steps".format(numpy.mean(sampler.get_autocorr_time(c=1,quiet=True))))
