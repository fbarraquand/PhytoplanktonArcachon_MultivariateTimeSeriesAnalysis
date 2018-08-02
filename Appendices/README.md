This folder contains all the scripts used in the Appendices (results and graphs). 

* Appendix 1 focuses on data analysis (time series plots to describe the environment at both our study sites, auto- and cross-correlation analyses in both time- and frequency-domain for both biotic and abiotic variables, interpolation for missing values).

* Appendix 2 focuses on both univariate and multivariate auto-regressive analyses, first by computing linear models with different lags on each time series, then by computing MAR(1) models. It also tests the ability of the MARSS package to retrieve the values of the parameters used in simulation (estimations are performed in the `Appendix2_analysis_simulated_data_scenario[1-3].r`, and corresponding diagnostics, i.e., number of false negatives, are drawn in `Appendix2_plot_for_simulated_data.r`.

* Appendix 5 contains searches for non-linearities in the model, such as thresholds in the environment or the community that would modify the dynamics of the two main groups of phytoplankton, or presence of a story effect.  
