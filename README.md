# PCDH19_FOOOF
 
This is the FOOOF algorithm developed by Haller et al. 2018, used in my PhD thesis on PCDH19 to find the slope of the power spectrum of LFP data.

All functions in this folder are developed by Haller et al., only the fooof_complete_didi.m script is my code, containing my implementation of it.
The program takes as an input .dat files containing the power spectrum of LFP recordings calculated with ZebraExplore (see my ZebraExplore repository). 
The program then calls FOOOF to fit an exponential function to the power spectra and find oscillations in the power spectra. It also provides measures 
of the goodness of fit of the exponential function. Finally, fooof_complete_didi.m calculates entropy as in Foffani et al. 2007. 

A detailed description of my implementation of FOOOF can be found in the comments on top of the fooof_complete_didi.m script. 