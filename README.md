# ww_paper

This repository has all code needed to recreate the analyses conducted in the paper [Semiparametric Inference of Effective Reproduction Number Dynamics from Wastewater Pathogen Surveillance Data](https://arxiv.org/abs/2308.15770). 
Models were fit in `Julia`, while simulation of synthetic data, and visualization of results was done in `R`. 
To set up the `Julia` project, follow the instructions [here](https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project).
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

## Setting up the `Julia` environment. 
The results from this project were generated using `Julia 1.8.5`, which can be downloaded [here](https://julialang.org/downloads/oldreleases/). 
If you want to use a more recent version of [`Julia`](https://julialang.org/downloads/), delete the [Manifest.toml](https://github.com/igoldsteinh/ww_paper/blob/main/Manifest.toml) file after you have cloned the repo. 
Once you have `Julia` installed, from the terminal, navigate to the project root directory then type `julia`. 
Your terminal will look like:
```
julia>
```
Now type `]`. Your terminal should now look like:
```
(@v1.8) pkg>
```
Then use the following commands
```
activate .
```
and 
```
instantiate
```
If you kept the `Manifest.toml`, the exact `Julia` package versions used to generate the paper results will be downloaded. 
If you did not, possibly newer versions of the packages will be downloaded instead. 
It would be surprising if newer versions of the packages led to different results. 
More information on `Julia` environments is available in the [Environments documentation](https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project).
When executing code from this repo, be mindful of your project directory; the `Julia` package versions are specific to the project environment. 
If you are executing code outside of this project, the packages you installed as part of the environment will not be available. 

## Quarto and Julia
The [vignettes folder](https://github.com/igoldsteinh/ww_paper/tree/main/vignettes) has two Quarto vignettes which condense the model fitting workflow into one [`Julia` vignette](https://github.com/igoldsteinh/ww_paper/blob/main/vignettes/fit_eirr_ww.qmd) that demonstrates how to fit the EIRR-ww model to the Los Angeles wastewater data via MCMC and one [`R` vignette](https://github.com/igoldsteinh/ww_paper/blob/main/vignettes/process_eirr_ww.qmd) that uses the results of the `Julia` vignette to visualize the saved MCMC results. 
We recommend starting with these vignettes, as they provide more detailed explanations of the code than the original scripts.

To execute the vignettes, we recommend using the IDE [VS Code](https://code.visualstudio.com) with the [`Julia`](https://code.visualstudio.com/docs/languages/julia) and [Quarto](https://quarto.org/docs/tools/vscode.html) extensions. 
Additional information on compiling `Julia` Quarto files is available [here](https://quarto.org/docs/computations/julia.html). 
You can also run each chunk of `Julia` code by copy pasting it into the REPL (the interactive Julia environment that opens when you type `julia` in the terminal).

## Model fitting workflow
The original workflow for the main models involves multiple files. As an example, to generate results from the the EIRR-ww model, use [fit_eirrc_closed.jl](https://github.com/igoldsteinh/ww_paper/blob/main/scripts/fit_models/fit_eirrc_closed.jl) to fit the model, then [eirrc_closed_generate_pp_and_gq.jl](https://github.com/igoldsteinh/ww_paper/blob/main/scripts/generate_quantities/eirrc_closed_generate_pp_and_gq.jl) to re-scale the posterior and generate posterior predictive values, finally [process_results_eirrc_closed.R](https://github.com/igoldsteinh/ww_paper/blob/main/scripts/process_results/process_results_eirrc_closed.R) creates tidy versions of the posterior and posterior predictive summaries.
When summarising results from multiple simulations, [summarise_eirrc_closed.R](https://github.com/igoldsteinh/ww_paper/blob/main/scripts/process_results/summarise_eirrc_closed.R) creates summary outputs. 
Similarly named files exist for all models used in the paper. 

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

## Fitting the Huisman model
To fit the Huisman model use [fit_huisman.R](https://github.com/igoldsteinh/ww_paper/blob/main/scripts/fit_models/fit_huisman.R). To visualize the results, use the code in [visualize_fit_to_LA_data.R](https://github.com/igoldsteinh/ww_paper/blob/main/scripts/visualize_results/visualize_fit_to_LA.R) and [visualize_frequentist_metrics.R](https://github.com/igoldsteinh/ww_paper/blob/main/scripts/visualize_results/visualize_frequentist_metrics.R). 

You'll need to install some packages, use the code below to do so:
```
install.packages(c("tidyverse",
                   "lubridate",
                   "patchwork",
                   "viridis",
                   "EpiEstim",
                   "zoo",
                   "tidybayes"))
```
## Fitting Epidemia and Rt-estim-gamma models
To fit the Epidemia and Rt-estim-gamma model, use [fit_epidemia.R](https://github.com/igoldsteinh/ww_paper/blob/main/scripts/fit_models/fit_epidemia.R) and [fit_estimgamma.R](https://github.com/igoldsteinh/ww_paper/blob/main/scripts/fit_models/fit_estimgamma.R) visualize the results using [visualize_fit_to_LA_data.R](https://github.com/igoldsteinh/ww_paper/blob/main/scripts/visualize_results/visualize_fit_to_LA.R).
Both models are written in `Stan`, installation instructions for `rstan` are available [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).
We also use the [`Epidema` package](https://imperialcollegelondon.github.io/epidemia/index.html), which is only available on Github. 

The following code installs `Epidemia`. 
```
#install.packages("devtools")
devtools::install_github("ImperialCollegeLondon/epidemia")
```
To install other needed packages, use the code below:
```
install.packages(c("brms",
                   "truncnorm",
                   "sdprisk",
                   "rstanarm"))
```
