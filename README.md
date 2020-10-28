# PCDH19_FOOOF

## General Information 
This is the FOOOF algorithm developed by Haller et al. 2018, and my adaptation of it.

## Use in PCDH19 project

Used in my PhD thesis on PCDH19 to find the slope of the power spectrum of LFP data.

## Usage
The program fooof_complete_didi.m is the main program to run, it calls all other functions in this folder when it needs them. 
It takes as an input .dat files containing the power spectrum of LFP recordings calculated with ZebraExplore (see my PCDH19_ZebraExplore repository). 
The program then calls FOOOF to fit an exponential function to the power spectra and find oscillations in the power spectra. It also provides measures 
of the goodness of fit of the exponential function. Finally, fooof_complete_didi.m calculates entropy as in Foffani et al. 2007. 

A detailed description of my implementation of FOOOF can be found in the comments on top of the fooof_complete_didi.m script. 

## Contributors
All functions in this folder are developed by Haller et al., only the fooof_complete_didi.m script is my code, containing my implementation of it.

## Dependencies
It requires a working installation of Python (later than 3.5), with numpy and scipy packages. 





