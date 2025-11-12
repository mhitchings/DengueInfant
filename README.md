# DengueInfant

Code, data, and model results accompanying "Immune-mediated protection and enhancement of dengue drives patterns of infant cases in Brazil"

RCode folder contains code to: 
1. Make figures contained in manuscripts describing spatiotemporal patterns in infant dengue in Brazil (make_data_figures_clean.R)
2. Fit models of infant dengue and severe dengue to surveillance data in Brazil (infantdengue_runstanmod_clean.R) using rstan files (fit_infantdengue_nolambdauncertainty_modeloption_clean.stan and fit_infantdengue_lambdauncertainty_modeloption_clean.stan)
3. Make plots of model fits and predictions (make_model_figures_clean.R)

stan_fits folder contains model parameters for different model assumptions

Note: Most models will take 8-16 hours to fit on a typical machine
