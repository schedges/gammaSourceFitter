# gammaSourceFitter
Given a background run, source simulation, and source run, fits simulation + bgnd to source to get conversion from arbitrary unit to keV and energy-dependent energy resolution parameters

## Overview
This code was originally designed to be used to calibrate liquid scintillator detectors. It is a python RooFit based code that requires a simulation, background run, and source run. A RooAddPdf is created combining the background run and simulation (modified by scaling and energy resolution parameters), and the result is fit to the source data. There are several ways of doing this, outlined below.

* **fitSim_floatAmplitude_keysPDF.py**: Steps through resolution parameters and conversion to ADC parameters in discrete steps, fitting only the amplitude of the simulation PDF. The NLL is computed for each step, and the lowest NLL value is returned. This version uses a RooKeysPDF to model the simulation data and background data. The advantage of this is that the fits are done in an unbinned fashion. The disadvantage is that there is potential bias introduced by the RooKeysPdf, so care should be taken to see that it does a good job describing the data it approximates. The interpolation parameters and mirroring options when generating the RooKeysPDF can be used to control this.
* **fitSim_floatAmplitude_fixedBinning.py**: Same as above, but a binned version instead of using a RooKeysPdf. Requires a pre-scaling factor to bring the simulation and data into approximately the same range. 
* **fitSim_floatAmplitudeAndConversion.py**: Steps through resolution parameters in discrete steps, for each one fitting the amplitude of the simulation PDF and the conversion parameters from an arbitrary unit to keV. There is an option to use either fixed binning (in which case a prescaling factor is needed, an estimate for the conversion from keV to arbitrary energy unit, that is within 50% of the true value), or variable binning (in which case appropriate binning needs to be selected over the range of the source spectrum and simulation spectrum such that the important details of both are captured).
* **fitSim_floatEverying.py**: Experimental, makes a RooProdPDF of a resolution function and simulation. Does a binned fit, requiring the input of a pre-scaling factor. 

There are a few other approaches that could be tried:
* Using RooMorphPDFs to determine scaling and resolution parameters
* Using fft convolution to introduce the scaling parameter
* Adding a new column with the transformed variable

## How to use the code:
The user must input the names of their TTrees, branches corresponding to channel number, energy, and time (the latter only applicable for data trees, not simulation trees). The code assumes the background data tree structure and source data tree structure follow the same structure. Ranges for resolution and energy conversion variables needs to be specified. Basically everything up to the "main code starts here" line needs to be adjusted to the specifics of your data. There are a few things that can be adjusted later in the code if necessary.

## Output of the code:
* A TTree containing all data about the fits including parameter values, errors, status, nll value, etc.
* A TTree contaiing the above for the fit with the lowest NLL values
* TCanvases of any fit that replaces the previous best fit.

## Notes/pitfalls:
* When specifying fit ranges and variable import ranges, care should be taken that the import range (lowerEnergyBound, upperEnergyBound) are greater than the fit range. Data outside the fit range needs to be imported, as it can be smeared into the fit range when the energy resolution parameters are applied.
* Thresholds are not modeled, so care should be taken to ensure the fit range is far enough away from threshold such that the fits are not being affected by threshold effects.
* Binning can affect fits, so studies should be done on how binning influences the fits and errors should be manually incorporated into the final parameter values. 
* RooKeysPDFs need to do a good job approximating the data over the fit range. This should be checked for background RooKeysPdf and first few simulation RooKeysPdfs. The Mirroring option, and interpolation number (last argument when generating the RooKeysPdfs) may need to be adjusted.
* By default warnings are reduced, it's probably a good idea to turn these back on in the RooFit settings section and the fitTo command.

## Known bugs:
* There is a known issues with root 6.18/00 and 6.18.02
* There is a memory leak where RooFit objects are not deleted in the loop. I've done what I can to minimize this, but think it is a RooFit bug.

## Theory:
### Liquid scintillators 
The light output of liquid scintillators is typically thought to be linear with energy at higher energies (>~50 keVee) [1]. Conversion from an arbitrary energy unit (peak high value, integral, max of an FIR, etc.) to a real energy unit (keV, MeV) is typically done using two parameters, a slope and an offset (reflecting low-energy non-linearities)

Light output = slope*(keV - offset)

where offsets range from ~0-25 keV. Typically in liquid scintillators, it is difficult to see full deposition peaks from gamma sources, so using Compton peaks to calibrate is often not possible. Fitting the Compton edge is possible, but energy resolution needs to be understood to do this reliably. Energy resolution is usually parameterized as

FWHM/E = 2.355 * sigma/E = sqrt(alpha^2 + beta^2/E + gamma^2)

 with typical values of alpha ranging from 0.1-0.15, beta  from 0.05-0.135, and gamma  from 0.002-0.065.

### Plastic Scintillators
The above typically holds true for plastic scintillators, but generally different energy resolution functions are used, either

sigma = 1/2.355 * sqrt(a * E)

with a between 0.014 and 0.032 [2], or

sigma = sqrt(A * E^2 + B * E).

with A around 0.0162, and B around 0.0105 [3]


### Sources
1. https://www.sciencedirect.com/science/article/pii/S0168900298006792
2. https://arxiv.org/abs/1004.3779
3. https://www.sciencedirect.com/science/article/pii/0029554X71903703
