README
================

## Running MCMC chains - scripts in MCMC folder

Test run of 3 chains x 1000 iterations each: - Run\_testchain1000.R

Production runs were conducted in 3 steps: (1) Run\_chain5000.R (first
5000 iterations) (2) Run\_chain2x5000.R (next 5000 iterations) (3)
Run\_chain4x5000.R (next 10000 iterations)

The results were summarized using the script: -
summarize\_mcmc\_results.R

## Updating posterior predictions with new testing data

The script update\_posterior\_predictions.R, located in MCMC folder,
downloads new testing data from COVID Tracking Project and generates new
posterior prediction
quantiles.

## Additional validation with international data - in International folder

The script get\_intl\_post\_pred.R downloads testing data from OWID and
UK governmen, and generates estimates of prevalence and seroprevalence
based on the posterior distributions from U.S. states.

## Manuscript figures - scripts in Figures folder

1)  Figure1AB\_Concept.R - Top 2 panels of Figure 1
2)  Figure2\_Seroprevalence.R - Figure 2
3)  Figure3\_Prevalence.R - Figure 3
4)  FIgure4\_map.R - Figure 4

## Supplemental Figures - scripts in Figures folder

1)  SupFig\_ParameterPosteriors.R
      - SupFig\_RandomEffects.pdf
2)  SupFig\_BayesianScatter.R
      - SupFig\_Bayesian\_scatter.pdf
3)  SupFig\_BayesianPrevalenceScatter.R
      - SupFig\_Bayesian\_prevalence\_scatter\_comb.pdf
4)  SupFig\_CurrentPrevSeroprev.R
      - SupFig\_CurrentEstimates.pdf
5)  SupFig\_BiasEstimates.R
      - SupFig\_Biases.pdf
6)  SupFig\_DivergencesStates.R
      - SupFigDivergences.pdf
7)  SupFig\_Intl\_post\_pred\_obs.R
      - SupFig\_Intl\_post\_RE\_nolog.pdf

## Quick summary

The R Markdown file Semi-Empirical-Prev\_Seroprev.Rmd provides a quick
summary of state-by-state and US overall results, compared to cumulative
cases and published seroprevalence data.

## Code for R shiny web app

The code for the R shiny web app is in the folder
COVID-19-Prevalence-and-Seroprevalence
