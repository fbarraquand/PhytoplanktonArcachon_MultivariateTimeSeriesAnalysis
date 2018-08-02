[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1326510.svg)](https://doi.org/10.5281/zenodo.1326510)


# Multivariate time series analyses of Arcachon phytoplankton - abiotic and biotic drivers of community dynamics
This repository contains the files and the data needed to reproduce the analyses in Barraquand et al. Oikos 2018 *Coastal phytoplankton community dynamics and coexistence driven by intragroup density-dependence, light and hydrodynamics*. 

### Data files 
* `B7_base.csv`, `Teychan_base.csv` are phytoplankton counts at Teychan and B7 (+ some additional variables)
* `meteo.Rdata`, `nao.Rdata`, `vent.Rdata` are respectively the meteorological, climate indices and wind datasets

### Functions
* `functions_global.r` defines a number of functions for data pre-processing
* `spectrum.r` and `spectrum_comparison.r` compute the Fourier spectra and coherence (correlation in the spectral domain), respectively
* `AR_test_order.r` searches for the order of the autoregressive model which best represents the univariate time series
* `eigenvalueB.r` checks the eigenvalues of the interaction matrix estimated with MAR(1) models
* `REPHY_aic_bic.r` prints the AICc and BIC corresponding to the different MAR(1) models previously estimated

### Model fitting
* `ARLM_function.r` fits the autoregressive linear models (AR(p) and ARX(p)) and plot Fig. 4 in the main text with the script `use_ARLM.r`
* `MARSS_plankton.r` fits the MAR(1) models on the plankton log(abundance) data
* `PlanktonLike_MAR_SimulationLoop.R` simulates plankton-looking data according to three scenarios (1:abiotic-forcing-only, 2:biotic-interactions-only, 3:both)
* `phase_models.r` fits the bivariate (plankton genus AST vs CHA), phase-dependent TAR models
* `Ricker.R` considers the alternative Ricker functional form of density-dependence (Appendix S5), and compare its fit and simulated dynamics to the Gompertz assumption (i.e., MAR(1) model in log(abundance)) 
* `arlm_aicc.csv` contains the different AICc for ARX(p), which are analysed in `figure_3.r` to produce Figure 3 in the main text
* The subrepo `MARSS_results` contains results from `MARSS_plankton.r` for B7 and Teychan for the five different interaction matrices considered
* `figure_5.r` and `figure_6.r` produce the corresponding figures in the main text

