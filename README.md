# gammaSourceFitter
Given a background run, source simulation, and source run, fits simulation + bgnd to source to get conversion from arbitrary unit to keV and energy-dependent energy resolution parameters

## Overview
This code was originally designed to be used to calibrate liquid scintillator detectors. It is a python RooFit based code that requires a simulation, background run, and source run. A RooAddPdf is created combining the background run and simulation (modified by scaling and energy resolution parameters), and the result is fit to the source data. There are several ways of doing this, outlined below.

* fitSim_floatAmplitude_keysPDF.py: Steps through resolution parameters and conversion to ADC parameters in discrete steps, fitting only the amplitude of the simulation PDF. The NLL is computed for each step, and the lowest NLL value is returned. This version uses a RooKeysPDF to model the simulation data and background data. The advantage of this is that the fits are done in an unbinned fashion. The disadvantage is that there is potential bias introduced by the RooKeysPdf, so care should be taken to see that it does a good job describing the data it approximates. The interpolation parameters and mirroring options when generating the RooKeysPDF can be used to control this.
* fitSim_floatAmplitude_fixedBinning.py: Same as above, but a binned version instead of using a RooKeysPdf. Requires a pre-scaling factor to bring the simulation and data into approximately the same range. 
* fitSim_floatAmplitudeAndConversion_variableBinning.py: Steps through resolution parameters in discrete steps, for each one fitting the amplitude of the simulation PDF and the conversion parameters from an arbitrary unit to keV. This version uses variable binning to allow capturing important details while keeping the binning from growing too large.
* fitSim_floatAmplitudeAndConversion_fixedBinning.py: Same as above, but requires a prescaling factor to bring the simulated data and measured data to approximately the same scale. This way the binning does not have to be extreme to capture the important details of the two data sets.
* fitSim_floatEverying.py: Experimental, makes a RooProdPDF of a resolution function and simulation. Does a binned fit, requiring the input of a pre-scaling factor. 

There are a few other approaches that could be tried:
* Using RooMorphPDFs to determine scaling and resolution parameters
* Using fft convolution to introduce the scaling parameter

## How to use the code:
The user must input the names of their TTrees, branches corresponding to channel number, energy, and time (the latter only applicable for data trees, not simulation trees). The code assumes the background data tree structure and source data tree structure follow the same structure.

## Output of the code:
* A ttree containing all data about the fits including parameter values, errors, status, nll value, etc.
* A ttree contaiing the above for the fit with the lowest NLL values
* TCanvases of any fit that replaces the previous best fit.

## Notes/pitfalls:
* When specifying fit ranges and variable import ranges, care should be taken that the import range (lowerEnergyBound, upperEnergyBound) are greater than the fit range. Data outside the fit range needs to be imported, as it can be smeared into the fit range when the energy resolution parameters are applied.
* Thresholds are not modeled, so care should be taken to ensure the fit range is far enough away from threshold such that the fits are not being affected by threshold effects.
* Binning can affect fits, so studies should be done on how binning influences the fits and errors should be manually incorporated into the final parameter values. 
* RooKeysPDFs need to do a good job approximating the data over the fit range. This should be checked for background RooKeysPdf and first few simulation RooKeysPdfs. The Mirroring option, and interpolation number (last argument when generating the RooKeysPdfs) may need to be adjusted.
