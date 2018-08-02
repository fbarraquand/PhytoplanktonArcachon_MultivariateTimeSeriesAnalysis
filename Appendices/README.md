This folder contains all the scripts used in the Appendices (results and graphs). 

* Appendix 1 focuses on data analysis (time series plots to describe the environment at both our study sites, auto- and cross-correlation analyses in time- and frequency-domain for both biotic and abiotic variables, interpolation for missing values).

* Appendix 2 focuses on both univariate and multivariate auto-regressive analyses. We first compute linear models with different lags on each time series, and then MAR(1) models. The code also tests the ability of the MAR(1) models to retrieve the values of the parameters used in simulation (estimations are performed in the `Appendix2_analysis_simulated_data_scenario[1-3].r`, and corresponding diagnostics, i.e., number of false negatives, are drawn in `Appendix2_plot_for_simulated_data.r`).

* Appendix 5 contains additional models with more non-linearities. In particular, we consider thresholds in the environmental variables or the community densities that would modify the dynamics of the two main groups of phytoplankton. 
