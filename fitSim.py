import ROOT
import sys
import numpy as np
import helperFns
import os
import emcee
import tqdm
import matplotlib.pyplot as plt
import corner
import csv as csvlib
import multiprocessing
import matplotlib

#Make sure process is single-threaded, we'll use multiprocessing library instead
os.environ["OMP_NUM_THREADS"] = "1"
nCores=multiprocessing.cpu_count()-1

# Reduces warning messages.
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
# Controls the number of steps used for integration. I needed to increase this
# in some fits when generating the NLL
ROOT.RooAbsReal.defaultIntegratorConfig().getConfigSection("RooIntegrator1D").setRealValue("maxSteps",30);
# Set integration step sizes, useful for balancing speed vs. accuracy
ROOT.RooAbsReal.defaultIntegratorConfig().setEpsAbs(1e-8)
ROOT.RooAbsReal.defaultIntegratorConfig().setEpsRel(1e-8)

#################
#Load input file#
#################
if not len(sys.argv)==2:
    print("Error! Wrong # of inputs")
    print("Usage: python fitSim.py <input file.dat>")
    sys.exit()

argDict=helperFns.parseFile(sys.argv[1])
fitRangeMin=argDict["fitRange"][0]
fitRangeMax=argDict["fitRange"][1]
binningType=argDict["binningType"]
nWalkers=argDict["nWalkers"]
nSteps=argDict["nSteps"]
nBurnInSteps=argDict["nBurnInSteps"]
title=argDict["title"]

###########################
#Create data & sim objects#
###########################
bgndObject=helperFns.dataObject(argDict["bgndName"],argDict["bgndFilePath"],argDict["bgndTreeName"],argDict["bgndEnergyBranchName"],argDict["bgndNormalization"],argDict["bgndNormalizationParameter"],argDict["bgndAdditionalCuts"],fitRangeMin,fitRangeMax)
#Sources
nSources=len(argDict["sourceNames"])
sourceObjects=[]
for i in range(0,nSources):
    sourceObjects.append(helperFns.dataObject(argDict["sourceNames"][i],argDict["sourcePaths"][i],
    argDict["sourceTreeName"][i],argDict["sourceEnergyBranchName"][i],argDict["sourceNormalization"][i],
    argDict["sourceNormalizationParameter"][i],argDict["sourceAdditionalCuts"][i],fitRangeMin,fitRangeMax))
#Sims
simObjects=[]
for i in range(0,nSources):
    simObjects.append(helperFns.simObject(argDict["simNames"][i],argDict["simPaths"][i],
    argDict["simTreeName"][i],argDict["simBranchName"][i],argDict["simBranchUnits"][i],
    argDict["maxSimEntries"],argDict["simAdditionalCuts"][i]))


########################
#Set-up blank histogram#
########################
#Make a histogram with passed in binning--we'll make copies of this for our plots, PDFs, etc.
if binningType=="fixed":
    hist=ROOT.TH1D("hist","",int(argDict["binningParameter"])-1,fitRangeMin,fitRangeMax)
else:
    hist=ROOT.TH1D("hist","",len(argDict["binningParameter"])-1,np.array(argDict["binningParameter"]))

########################################
#Are we fitting or guessing parameters?#
########################################
#See if we want the user to guess at parameters (using first source/sim):
if argDict["runType"]=="guess":
    #Scale source hist by bgnd normalization
    sourceHist=hist.Clone("sourceHist")
    sourceHist=helperFns.fillHist(sourceHist,sourceObjects[0].data)
    sourceHist.Scale(bgndObject.normalization/sourceObjects[0].normalization)
    
    #Subtract off bgnd
    bgndHist=hist.Clone("bgndHist")    
    bgndHist=helperFns.fillHist(bgndHist,bgndObject.data)
    sourceHist.Add(bgndHist,-1)
    
    #Launch guessing program
    helperFns.guessParams(simObjects[0].data,sourceHist)

######################################
#See what's fixed and what's floating#
######################################
parNames=["alpha","beta","gamma","slope","offset"]
parameters=[] #Hold's a dictionary for each parameter
ndim=0
for parName in parNames:
    parameters.append(dict())
    parameters[-1]["name"]=parName
    if len(argDict[parName+"Range"])==2:
        print(parName+" is a free parameter")
        parameters[-1]["status"]="free"
        parameters[-1]["index"]=ndim
        parameters[-1]["min"]=argDict[parName+"Range"][0]
        parameters[-1]["max"]=argDict[parName+"Range"][1]
        ndim+=1
    else:
        print(parName+" is a fixed parameter")
        parameters[-1]["status"]="fixed"
        parameters[-1]["value"]=argDict[parName+"Range"][0]

#Add one dimension for each source to represent amplitude
for sourceNum in range(0,nSources): 
    parameters.append(dict())
    parameters[-1]["name"]=argDict["sourceNames"][sourceNum]+"Scaling"
    if len(argDict["simScalingRange"])==2:
        print(argDict["sourceNames"][sourceNum]+"Scaling"+" is a free parameter")
        parameters[-1]["status"]="free"
        parameters[-1]["index"]=ndim
        parameters[-1]["min"]=argDict["simScalingRange"][0]
        parameters[-1]["max"]=argDict["simScalingRange"][1]
        ndim+=1
    else:
        print(argDict["sourceNames"][sourceNum]+"Scaling"+" is a fixed parameter")
        parameters[-1]["status"]="fixed"
        parameters[-1]["value"]=argDict["simScalingRange"][0]
         
###############
#RooFit Set-up#
###############
energyVar = ROOT.RooRealVar("energyVar","energyVar",fitRangeMin,fitRangeMax)
# Make RooArgSet, can lead to memory leaks if made in a loop
argSet=ROOT.RooArgSet(energyVar)
argList=ROOT.RooArgList(energyVar)

#Make binning
binning=ROOT.RooBinning(fitRangeMin,fitRangeMax)
if binningType=="fixed":
    bins=np.linspace(fitRangeMin,fitRangeMax,int(argDict["binningParameter"]))
else:
    bins=[i for i in argDict["binningParameter"]]
for i in bins:
    binning.addBoundary(i)
energyVar.setBinning(binning)

#Load source data object to a RooDataSet
dataSets=[]
for i in range(0,nSources):
    dataSets.append(ROOT.RooDataSet.from_numpy({"energyVar":sourceObjects[i].data}, [energyVar]))
#Load bgnd data into a RooHistPdf
bgndDataSet = ROOT.RooDataSet.from_numpy({"energyVar":bgndObject.data}, [energyVar])
bgndDataHist = bgndDataSet.binnedClone()
bgndHistPdf = ROOT.RooHistPdf("bgndHistPdf","bgndHistPdf",argSet,bgndDataHist,0) #1 specifies interpolation order

#Calculate expected # of bgnd counts in each run
expectedBackgroundCounts=[]
for i in range(0,nSources):
    expectedBackgroundCounts.append(bgndObject.countsInFitRange*sourceObjects[i].normalization/float(bgndObject.normalization))
#TODO: dead time calculator. 
#source run rate - bgnd run rate = increased source rate. 
#Increased source rate * waveform length = fractional dead time
#Expected bgnd rates (w/dead time) = expectedBackgroundCounts*(1-fractional dead time)
   
expectedSourceCounts=[]
for i in range(0,nSources):
    expectedSourceCounts.append(sourceObjects[i].countsInFitRange-expectedBackgroundCounts[i])
    print("Expected "+str(expectedSourceCounts[-1])+" source counts")
################
##EMCEE Set-up##
################
PLOT=0

#Make into numpy arrays 
labels=[i["name"] for i in parameters if i["status"]=="free"]
mins=[i["min"] for i in parameters if i["status"]=="free"]
maxes=[i["max"] for i in parameters if i["status"]=="free"]
pos_min = np.array(mins)
pos_max = np.array(maxes)

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
pos = [pos_min + psize*np.random.rand(ndim) for i in range(nWalkers)]

######################EMCEE FUNCTIONS DEFINED HERE##############################

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
    return -np.inf
  
#This is our RooFit function that returns the POSITIVE log likelihood
def lnlike(theta):

  # Reset nllVal to zero--we add NLLs from each source fit
  nllVal=0
  
  #Load theta values into parameters
  alpha=theta[parameters[0]["index"]] if parameters[0]["status"]=="free" else parameters[0]["value"]
  beta=theta[parameters[1]["index"]] if parameters[1]["status"]=="free" else parameters[1]["value"]
  gamma=theta[parameters[2]["index"]] if parameters[2]["status"]=="free" else parameters[2]["value"]
  slope=theta[parameters[3]["index"]] if parameters[3]["status"]=="free" else parameters[3]["value"]
  offset=theta[parameters[4]["index"]] if parameters[4]["status"]=="free" else parameters[4]["value"]
  
  for i in range(0,nSources):
    #Make ampltitude vars
    sourceScaling=theta[parameters[i+5]["index"]] if parameters[i+5]["status"]=="free" else parameters[i+5]["value"]
    sourceCountsVar = ROOT.RooRealVar("sourceCountsVar","sourceCountsVar",sourceScaling*expectedSourceCounts[i])
    sourceCountsVar.setConstant(1)
    
    bgndCountsVar = ROOT.RooRealVar("bgndCountsVar","bgndCountsVar",expectedBackgroundCounts[i])
    bgndCountsVar.setConstant(1)
    
    #Make smeared sim PDF 
    smearedSimArray=helperFns.smearAndShiftArray(simObjects[i].data,alpha,beta,gamma,slope,offset)
    tempHist=hist.Clone("tempHist")
    tempHist=helperFns.fillHist(tempHist,smearedSimArray)
    smearedSimDataHist=ROOT.RooDataHist("smearedSimDataHist","smearedSimDataHist",argList,tempHist)
    simPdf = ROOT.RooHistPdf("simPdf","simPdf",argSet,smearedSimDataHist,0) #1 specifies interpolation order
    
    ##Make Model
    pdfList = ROOT.RooArgList(bgndHistPdf,simPdf)
    ampList = ROOT.RooArgList(bgndCountsVar,sourceCountsVar)
    model = ROOT.RooAddPdf("model","model",pdfList,ampList)

    #Compute nll
    nll = model.createNLL(dataSets[i],
      ROOT.RooFit.Extended(1),
      ROOT.RooFit.Verbose(0),
      ROOT.RooFit.NumCPU(1)
    )
    ROOT.SetOwnership(nll,True)
    
    #Make NLL positive
    nllVal += (-1*nll.getVal())

      
    if PLOT==1:     
        modelHist=hist.Clone("modelHist")
    
        sourceDataHist=dataSets[i].binnedClone()
        sourceHist=ROOT.RooAbsData.createHistogram(sourceDataHist,"energyVar",energyVar,ROOT.RooFit.Binning("binning"))
        sourceCounts=sourceObjects[i].countsInFitRange
        
        bgndHist=ROOT.RooAbsData.createHistogram(bgndDataHist,"energyVar",energyVar,ROOT.RooFit.Binning("binning"))
        bgndCounts=expectedBackgroundCounts[i]
        
        simCounts=sourceScaling*expectedSourceCounts[i]
        helperFns.plotResults(sourceHist,sourceCounts,bgndHist,bgndCounts,tempHist,simCounts,modelHist,title+"_fit_"+sourceObjects[i].name)
        del sourceDataHist
        del sourceHist
        del bgndHist
      
    #Memory management
    #Assume fixed?
    nll.Delete()
    del nll
    del model
    del ampList
    del pdfList
    del simPdf
    del smearedSimDataHist
    del tempHist
    del bgndCountsVar
    del sourceCountsVar
    
  #Return total nllval from all sources
  return nllVal
  
#Calls out roofit function, makes sure output is not infinite and parameters in allowed range 
def lnprob(theta):
  lp = lnprior(theta)
  if not np.isfinite(lp):
    return -np.inf
  return lp + lnlike(theta)
#emcee call


##########################MC RUN STARTS HERE####################################
with multiprocessing.Pool(nCores) as pool:
  sampler = emcee.EnsembleSampler(nWalkers,ndim,lnprob,pool=pool)
  #Burn in
  print("Starting burn in...")
  pos, prob, state  = sampler.run_mcmc(pos, nBurnInSteps, progress=True)
  print("Burn-in complete! Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
  sampler.reset()
  pos, prob, state  = sampler.run_mcmc(pos,nSteps,progress=True)


###########################MC PLOTTING HERE#####################################
#Save plots as a PDF
matplotlib.use('PDF')

#Get the LL values
samples=sampler.flatchain
lnprobs = sampler.lnprobability[:,:]
flatLnprobs = lnprobs.reshape(-1)
with open("sampler_"+title+".csv", 'w') as samplerOutFile:
  theWriter = csvlib.writer(samplerOutFile, delimiter=',')
  for sampleLine, llval in zip(samples, flatLnprobs):
    theWriter.writerow(np.append(sampleLine,llval))

#MAKE TRACE PLOTS OF EACH PARAMATER
fig = plt.figure(figsize=(10,ndim*2))
gs = fig.add_gridspec(ndim,1)
plt.subplots_adjust(hspace=0.4)
for i in range(0,ndim):
  axes = fig.add_subplot(gs[i,:])
  axes.plot(sampler.chain[:,:,i].T, '-', color='k', alpha=0.3)
  axes.set_title(labels[i])
plt.savefig("tracePlot_"+title+".pdf")

#CORNER PLOT HERE 
samples=sampler.flatchain
fig = corner.corner(samples, labels=labels, ranges=ranges, quantiles=[0.16,0.5,0.84],show_titles=True,title_kwargs={'fontsize':12})
fig.savefig("corner_"+title+".pdf")

#CALCULATE QUANTILES HERE
bestFitValues=[]
for i in range(ndim):
  mcmc = np.percentile(sampler.chain[:,:,i],[16, 50, 84])
  q = np.diff(mcmc)
  print(labels[i]+": "+str(mcmc[1])+"+"+str(q[0])+" -"+str(q[1])+"\n")
  bestFitValues.append(mcmc[1])

#Plot best fit 
PLOT=1
lnlike(bestFitValues)


#Print out stats
print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction))+"\n")
print("Mean autocorrelation time: {0:.3f} steps".format(np.mean(sampler.get_autocorr_time(c=1,quiet=True))))
#Plot and save results
