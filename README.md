# Multivariate time series analyses of Arcachon phytoplankton time series - abiotic and biotic drivers of community dynamics
This repository contains the files and the data needed to reproduce the analyses in Barraquand et al. 2017 *Weak interactions between groups and predominantly physical drivers of community dynamics in coastal phytoplankton*. Specifically, 

* `B7_base.csv`, `Teychan_base.csv` are phytoplankton counts at Teychan and B7 (+ some additional variables)
* `meteo.Rdata`, `nao.Rdata`, `vent.Rdata` are respectively the meteorological, climate indices and wind datasets
* `functions_global.r` defines a number of functions for data pre-processing
* `spectrum.r` and `spectrum_comparison.r` compute the Fourier spectra and coherence (correlation in the spectral domain), respectively
* `ARLM_fixed_nblag_clean.r` fits the autoregressive linear models (AR(p) and ARX(p)) and plot Fig. 4 in the paper
* `MARSS_plankton.r` fits the MAR(1) models on the plankton log(abundance) data
* `PlanktonLike_MAR_SimulationLoop.R` simulates plankton-looking data according to three scenarios (1:abiotic-forcing-only, 2:biotic-interactions-only, 3:both)
* `phase_models.r` fits the bivariate (plankton genus AST vs CHA), phase-dependent TAR models
* `Ricker.R` considers the alternative Ricker functional form of density-dependence, and compare its fit and simulated dynamics to the Gompertz assumption (i.e., MAR(1) model in log(abundance)) 
* `arlm_aicc.csv` contains the different AICc for ARX(p), which are analysed in `figure_3.r` to produce Figure 3 in Barraquand et al. 2017
* `MARSS_results` contains results from `MARSS_plankton.r` for B7 and Teychan for the five different interaction matrices considered
* `figure_5.r` and `figure_6.r` produce the corresponding figures in the main text

