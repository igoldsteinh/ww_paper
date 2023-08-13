# ww_paper

This repository has all code needed to recreate the analyses conducted in the paper Semiparametric Inference of Effective Reproduction Number Dynamics from Wastewater Pathogen Surveillance Data. 
Models were fit in Julia, while simulation and visualization was done in R. 
[R](https://github.com/igoldsteinh/testpackage) and [Julia](https://github.com/igoldsteinh/testpackage.jl) packages which implement the models described in the paper are also available. 

## Navigation
Core functions--including simulation engines and implementations of Bayesian models--are located in the [src folder](https://github.com/igoldsteinh/ww_paper/tree/main/src). 
Scripts for simulating data, fitting models, processing results, and visualizing results, are located in the [scripts folder](https://github.com/igoldsteinh/ww_paper/tree/main/scripts). 
As an example, to fit the EIRR-ww model to the data from Los Angeles, CA, set `sim=="real"` for each of the next three files and execute [scripts/fit_models/fit_eirrc_closed.jl](https://github.com/igoldsteinh/ww_paper/blob/main/scripts/fit_models/fit_eirrc_closed.jl) to fit the model. 
Then [scripts/generate_quantities/eirrc_closed_generate_pp_and_gq.jl](https://github.com/igoldsteinh/ww_paper/blob/main/scripts/generate_quantities/eirrc_closed_generate_pp_and_gq.jl) to re-scale the posterior and generate posterior predictive values. 
Then use [scripts/process_results/process_results_eirrc_closed.R](https://github.com/igoldsteinh/ww_paper/blob/main/scripts/process_results/process_results_eirrc_closed.R) to create tidy versions of the posterior and posterior predictive summaries.
The workflow is similar for fitting to simulated data (simply change the `sim` and `seed` value), or for using a different model (use the correspondingly named file).
