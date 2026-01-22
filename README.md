# R code for Bachelor AM Thesis simulations

This repository contains the R code used for the simulation studies in my Bachelor's thesis in Applied Mathematics at the University of Twente.
The code corresponds with the method described in the thesis and is made publicly available to ensure transparency and reproducibility.

## Overview of files
- "Ettest.R"
  Implements the e-value t-test functions 'e.t.test1' and 'e.t.test2', which are used
  throughout the e-value simulation studies.

- "e_value_simulation.R"
  Runs the main simulation study for the sequential e-value t-test.
  For a range of true standardised effect sizes, the script estimates power and stopping sample sizes for two e-value variants:
  - e1: estimated effect size
  - e2: oracle effect size  
  The script also produces the figures used in the thesis and poster.

- "powersimulations.R"
  Determines the maximum sample size ('nmax') required to achieve the target power for the e-value procedure, based on repeated sequential simulations.

- "SD_variations.R"
  Contains the pilot-plus-t-test simulation study.
  Different methods for estimating or stabilising the standard deviation in the pilot phase are compared, and their impact on planned sample sizes is investigated.

- "geen12.R"
  Investigates the instability of the estimated standardised effect size when the pilot sample size is small.
  This script illustrates the sampling variability of the effect size estimator and supports the discussion on the unreliability of pilot-based planning.

## Reproducibility
All simulation scripts set a fixed random seed to ensure reproducibility.

## Notes
Some scripts assume the presence of auxiliary input files (e.g. tables with 'nmax' values) or source functions defined in other scripts. 
Paths may need to be adjusted depending on the local setup.
