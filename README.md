# ww_paper

This repository has all code needed to recreate the analyses conducted in the paper Semiparametric Inference of Effective Reproduction Number Dynamics from Wastewater Pathogen Surveillance Data. 
Models were fit in `Julia`, while simulation of synthetic data, and visualization of results was done in `R`. 
To set up the `Julia` project, follow the instructions [here](https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project).
If you want to execute this code using a different version of `Julia` than originally used, delete the [Manifest.toml](https://github.com/igoldsteinh/ww_paper/blob/main/Manifest.toml) file.
[`R`](https://github.com/igoldsteinh/concRt) and [`Julia`](https://github.com/igoldsteinh/concRt.jl) packages which implement the models described in the paper are also available. 
All results files needed to reproduce the figures are in the repo, individual simulation results are excluded for the sake of storage.

## Navigation
```
├── data                          <- Processed real and simulated data
│   └── sim_data                  <- Simulated data
│
├── figures                       <- Paper figures
│
├── raw_data                      <- Unprocessed real data
│
├── results                       <- Model outputs, organized by model, and then by output type
│   │                                Example for the EIRR-ww model (eirrc) is shown
│   │                                Structure is the same for main models
│   └── eirrc                     <- Summaries of simulation results
│        ├── generated_quantities <- Correctly scaled posteriors, quantiles and mcmc samples
│        ├── posterior_predictive <- Posterior predictive mcmc samples and quantiles
│        └── posterior_samples    <- Raw Julia posterior samples
│
├── scripts                       <- Paper code 
│   ├── fit_models                <- Fit models to real and simulated data
│   ├── fit_splines               <- Fit splines to data for prior elicitation
│   ├── generated_quantities      <- Turn raw Julia MCMC output into more useable csv files
│   ├── process_real_data         <- Turn raw real data into processed real data
│   ├── process_results           <- Final processing of mcmc/summarising simulation results
│   ├── simulated_data            <- Simulate data/test simulation engines against each other
│   └── visualize_results         <- Turn summaries of model results into paper figures
├── slurm_submissions             <- Slurm scheduler files for use on computing cluster
│   
├── src                           <- Models, priors, simulation engines, utility functions
│   
├── vignettes                     <- Example code for fitting models and processing results
└──     
```

## Model fitting workflow
The workflow for the main models involves multiple files. As an example, to generate results from the the EIRR-ww model, [fit_eirrc_closed.jl](https://github.com/igoldsteinh/ww_paper/blob/main/scripts/fit_models/fit_eirrc_closed.jl) is used to fit the model, then [eirrc_closed_generate_pp_and_gq.jl](https://github.com/igoldsteinh/ww_paper/blob/main/scripts/generate_quantities/eirrc_closed_generate_pp_and_gq.jl) to re-scale the posterior and generate posterior predictive values, finally [process_results_eirrc_closed.R](https://github.com/igoldsteinh/ww_paper/blob/main/scripts/process_results/process_results_eirrc_closed.R) creates tidy versions of the posterior and posterior predictive summaries.
When summarising results from multiple simulations, [summarise_eirrc_closed.R](https://github.com/igoldsteinh/ww_paper/blob/main/scripts/process_results/summarise_eirrc_closed.R) creates summary outputs. 
Similarly named files exist for all models used in the paper. 
The [vignettes folder](https://github.com/igoldsteinh/ww_paper/tree/main/vignettes) has two vignettes which condense this workflow and demonstrate how to fit the EIRR-ww model to the Los Angeles wastewater data.  

## Simulation name key
When executing scripts, the `sim` parameter controls what simulation is being used, the `seed` parameter controls the seed and also the specific data set used. 
For simulations, we used values of `seed` from 1 to 100. 
For the analysis of the Los Angeles wastewater data, `seed=1`. 
Here is a key translating the values of `sim`:
* `sim=1` = `Baseline`
* `sim=3` = `10-rep`
* `sim=4` = `3-mean`
* `sim=5` = `10-mean`
* `sim=6` = `1-rep`
* `sim=8` = `Low Prop`
* `sim=9` = `Low Init`
* `sim=10` = `High Init`
* `sim="real"` = `Los Angeles JWPCP wastewater data`
* `sim="ODE"` = `Baseline data observed every 12 hours`

We do not speak of simulations 2 and 7. 

## Model name key
The model names used in the code are not the same as those used in the paper. 
Here is a key:
* `eirrc_closed` = `EIRR-ww`
* `eir_cases` = `EIR-cases`
* `seir_cases` = `SEIR-cases`
* `seirr_student` = `SEIRR-ww`
* `huisman` = `Huisman`
* `epidemia` = `Epidemia`
* `estimgamma` = `Rt-estim-gamm`
* `eirr` = `EIRR-ww with ODE solver` 
