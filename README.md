README
================

## Folders

### Root
- Readme file
- setup_data.R script (contains switches for different model formulations)

### Data 
- Testing data 
- Seroprevalence survey data
- Basic data on states (population, FIPS, long/lat for map, etc.)

### Figures
- Scripts for making figures
- Figure PDF and JPEG files

### Functions
- Scripts defining functions
- Scripts running analyses

### International
- Scripts and data for doing international comparisons

### MCMC
- Results of primary model fitting

### MCMC.n0.X
- Results of model fitting for fixed n=0.X

### MCMC.tauX
- Results of model fitting for alternative tau=X

## Running MCMC chains - scripts in MCMC folder

First modify the "switches" in the "Model options" section of setup_data.R to conform to the model formulation you want to run.

Shorter runs (e.g., for fixed n) are performed by running script Run_4chains1000x2.R. This will run 4 chains of length 1000, then restart those chains to run for another 1000.  Then, the first 200 iterations are discarded, and the results summarized by the script summarize_mcmc_results.R.

Longe rruns (e.g., for random effects n) are performed by running each chain separately using the scripts Run_onechain4x5000.X.R, where X=1...4.  For each chain, this will run 5000 iterations, restart for another 5000, then restart again for 10000 iterations.  After all 4 chains are completed, the summary script summarize_mcmc_results.R needs to be run "manually."

The script summarize_mcmc_results.R generates random subsample for analysis, and then calculates new posterior prediction quantiles.

## Additional validation with international data - in International folder

The script get_intl_post_pred.R uses testing data from OWID and UK government, and generates estimates of prevalence and seroprevalence based on the posterior distributions from U.S. states.

## Manuscript figures - scripts in Figures folder

1. Figure1BC_Concept.R - Lower 2 panels of Figure 1
    - Figure1BC_Concept.pdf/jpeg
2. Figure2_Seroprevalence.R - Figure 2
    - Figure2_Seroprevalence_posteriors_validation_page.pdf/jpeg
3. Figure3_Prevalence.R - Figure 3
    - Figure3_v1_prevalence_posteriors_page.pdf/jpeg
4. Figure4_map.R - Figure 4
    - Figure4_Maps.pdf/jpeg
    
## Supplemental Figures - scripts in Figures folder ()

- (S1) SupFig_ParameterPosteriors.R 
    - SupFig_RandomEffects.pdf/jpeg 
- (S2) SupFig_BayesianScatter.R
    - SupFig_Bayesian_scatter.pdf/jpeg 
- (S3) SupFig_BayesianScatter_Validation.R
    - SupFig_Bayesian_scatter_validation.pdf/jpeg 
- (S4) SupFig_BayesianPrevalenceScatter.R
    - SupFig_Bayesian_prevalence_scatter_comb.pdf/jpeg 
- (S5) SupFig_CurrentPrevSeroprev.R
    - SupFig_CurrentEstimates.pdf/jpeg 
- (S6) Created with Figure4_map.R
    - Figure4_Maps.MCMC.n0.5.pdf/jpeg
- (S7) Created with SupFig_CurrentPrevSeroprev.R
    - SupFig_CurrentEstimates.MCMC.n0.5
- (S8) SupFig_BiasEstimates.R
    - SupFig_Biases.pdf/jpeg 
- (S9) SupFig_DivergencesStates.R
    - SupFigDivergences.pdf/jpeg 
- (S10) SupFig_Intl_post_pred_obs.R
    - SupFig_Intl_post_RE_nolog.pdf/jpeg 
- (S11-S14) SupFig_TauSensitivity.R
    - SupFig_TauSensitivity_parameters.pdf/jpeg
    - SupFig_TauSensitivity_SP.pdf/jpeg
    - SupFig_TauSensitivity_I.pdf/jpeg
    - SupFig_TauSensitivity_Itot.pdf/jpeg

## Quick summary

The R Markdown file Semi-Empirical-Prev_Seroprev.Rmd provides a quick summary of state-by-state and US overall results, compared to cumulative cases and published seroprevalence data.

## Code for R shiny web app

The code for the R shiny web app (https://wchiu.shinyapps.io/COVID-19-Prevalence-and-Seroprevalence/) is in the folder 

- COVID-19-Prevalence-and-Seroprevalence
