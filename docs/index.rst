Gamma Source Fitter Documentation
=================================
https://github.com/schedges/gammaSourceFitter

Documentation for python code to fit gamma calibration sources given source data, background data, and a sim. Includes required packages, basic usage, physics overview, and troubleshooting advice.

.. toctree::
   :maxdepth: 2
   
   physics
   requirements
   usage
   troubleshooting
   support
   
## Overview
This code was originally designed to be used to calibrate liquid scintillator detectors using gamma sources. It is a python and RooFit-based code that requires a simulation of the source (MCNP, GEANT, etc.), a background run without source, and a source run. A RooAddPdf is created combining the background run and simulation (modified by scaling and energy resolution parameters that float), and the result is fit to the source data.

* **fitSim_emcee_simultaneous.py**: Uses the python emcee package (https://emcee.readthedocs.io/en/stable/). The recommended way to use this code. Requires matplotlib, tqdm, emcee, and corner. Can fit multiple gamma sources simultaneously.
* **fitSim_floatAmplitude.py**: Non-emcee version of the code, not recommended to use! Steps through resolution parameters and energy conversion parameters in discrete steps, fitting only the amplitude of the simulation PDF for each set of parameters. The NLL is computed for each step, and the lowest NLL value is returned. Can only fit a single gamma source at a time.

## How to use the code:
 The user must specify some set-up at the start of the code (binning, background pdf type, file paths, tree/branch names, channel numbers, etc.) The code assumes the background data tree structure and source data tree structure follow the same structure. Ranges for resolution and energy conversion variables needs to be specified. 

## Output of fitSim_emcee_simultaneous.py:
* Sampler output (csv file)
* Corner plots with quantiles
* Trace plots
* Plot of best fit parameters

## Output of fitSim_floatAmplitude.py:
* A TTree containing all data about the fits including parameter values, status, nll value.
* Profile likelihood TGraph of each parameter, used for getting errors.
* Plot of the background + simulation + source data for the best fit parameters.

## Notes/pitfalls:
* When specifying fit ranges and variable import ranges, care should be taken that the import range (lowerEnergyBound, upperEnergyBound) are greater than the fit range. Data outside the fit range needs to be imported, as it can be smeared into the fit range when the energy resolution parameters are applied.
* Thresholds are not modeled, so care should be taken to ensure the fit range is far enough away from threshold such that the fits are not being affected by threshold effects.
* Binning and fit range are nuisance parameters, so some studies should be done on how these quantities affect fit results and errors. 
* If the background PDF type is a RooKeysPDF, care should be taken to check this does a good job approximating the background data. The Mirroring option, and interpolation number (last argument when generating the RooKeysPdfs) may need to be adjusted.
* If background PDF type is "binned", then you will get errors if you have too many bins such that the generated PDF is zero for some values
* It can be difficult to fit an offset with few or simple gamma sources used, so you may have to fix it in some instances. See theory notes below for more info.

## Known issues:
* There is a known issues with root 6.18/00 and 6.18.02
* There is an issue with RooFit when the createNll function is called too many times, resulting in a slowdown. You can use the emcee h5py backend to get around this by breaking up your call into multiple ones.. https://root-forum.cern.ch/t/non-linear-increase-of-execution-time-when-generating-more-toys/27411

## Theory:
### Liquid scintillators 
The light output of liquid scintillators is typically thought to be linear with energy at higher energies (>~50 keVee) [1]. Conversion from an arbitrary energy unit (peak high value, integral, max of an FIR, etc.) to a real energy unit (keV, MeV) is typically done using two parameters, a slope and an offset (reflecting low-energy non-linearities)

Light output [ADC] = slope [ADC/keV] * (energy [keV] - offset [keV])

where the offset is a small positive value, typically ~5 keV ([2],[3]). Typically in liquid scintillators, it is difficult to see full deposition peaks from gamma sources, so using Compton peaks to calibrate is often not possible. Fitting the Compton edge is possible, but energy resolution needs to be understood to do this reliably. Energy resolution is usually parameterized as

FWHM/E = 2.355 * sigma/E = sqrt(alpha^2 + beta^2/E + gamma^2)

 with typical values of alpha ranging from 0.1-0.15, beta  from 0.05-0.135, and gamma  from 0.002-0.065.

### Plastic Scintillators
The above typically holds true for plastic scintillators, but generally different energy resolution functions are used, either

sigma = 1/2.355 * sqrt(a * E)

with a between 0.014 and 0.032 [4], or

sigma = sqrt(A * E^2 + B * E).

with A around 0.0162, and B around 0.0105 [5]


### Sources
1. https://www.sciencedirect.com/science/article/pii/S0168900298006792
2. https://www.sciencedirect.com/science/article/pii/S0168900201014103
3. https://www.sciencedirect.com/science/article/pii/016890029290824N
4. https://arxiv.org/abs/1004.3779
5. https://www.sciencedirect.com/science/article/pii/0029554X71903703

