# Repo for manuscript "A quantum parallel Markov chain Monte Carlo"

## Code/

`analyzeBlackHole.R`: takes outputs from `blackHoleExample.R` and creates manuscript Figure 8

`analyzeIsing.R`: takes outputs from `IsingExp.R` and creates Figure 7

`blackHoleExample.R` runs Bayesian image analysis for black hole intensity map application

`groverConvergence.R` contains code for Grover's algorithm convergence experiment, creates Figure 1

`IsingExp.R` runs MCMC for 500-by-500 Ising target for number of proposals specified via command line

`multipleExperiments.R` contains code for creating Figures 2-6 and Table 1; to create the latter, code reads in results from `qpMCMC.ipynb`

`qpMCMC.ipynb` implements experiment on massive mixture of Gaussians, output of which contributes to Table 1

`runIsings.sh` bash script for running `IsingExp.R` for all proposal counts considered
