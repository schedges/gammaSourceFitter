## Overview
This code was originally designed to be used to calibrate liquid scintillator detectors using gamma sources. It is a python and RooFit-based code that requires a simulation of the source (MCNP, GEANT, etc.), a background run without source, and a source run. A RooAddPdf is created combining the background run and simulation (modified by scaling and energy resolution parameters that float), and the result is fit to the source data.

* **fitSim_emcee_simultaneous.py**: Uses the python emcee package (https://emcee.readthedocs.io/en/stable/). The recommended way to use this code.
* **fitSim_floatAmplitude.py**: Non-emcee version of the code, outdated! Steps through resolution parameters and energy conversion parameters in discrete steps, fitting only the amplitude of the simulation PDF for each set of parameters. The NLL is computed for each step, and the lowest NLL value is returned.

## How to use the code:
 The user must specify some set-up at the start of the code (binning, background pdf type, file paths, tree/branch names, channel numbers, etc.) The code assumes the background data tree structure and source data tree structure follow the same structure. Ranges for resolution and energy conversion variables needs to be specified. 

## Output of fitSim_floatAmplitude.py:
* A TTree containing all data about the fits including parameter values, status, nll value.
* Profile likelihood TGraph of each parameter, used for getting errors.
* Plot of the background + simulation + source data for the best fit parameters.

## Output of fitSim_emcee_simultaneous.py:
* Sampler output (csv file)
* Corner plots with quantiles
* Trace plots
* Plot of best fit parameters

## Notes/pitfalls:
* When specifying fit ranges and variable import ranges, care should be taken that the import range (lowerEnergyBound, upperEnergyBound) are greater than the fit range. Data outside the fit range needs to be imported, as it can be smeared into the fit range when the energy resolution parameters are applied.
* Thresholds are not modeled, so care should be taken to ensure the fit range is far enough away from threshold such that the fits are not being affected by threshold effects.
* Binning and fit range are nuisance parameters, so some studies should be done on how these quantities affect fit results and errors. 
* If the background PDF type is a RooKeysPDF, care should be taken to check this does a good job approximating the background data. The Mirroring option, and interpolation number (last argument when generating the RooKeysPdfs) may need to be adjusted.

## Known bugs:
* There is a known issues with root 6.18/00 and 6.18.02

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
