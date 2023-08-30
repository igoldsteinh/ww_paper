using DrWatson
using Revise
using JLD2
using FileIO
using CSV
using DataFrames
using Turing
using DifferentialEquations
using LogExpFunctions
using Random
using ForwardDiff
using Optim
using Random
using LineSearches
using ww_paper


sim =
if length(ARGS) == 0
  1
else
  parse(Int64, ARGS[1])
end

priors_only = sim == 0

seed = 
if length(ARGS) == 0
  1
else 
  parse(Int64, ARGS[2])
end 


if priors_only
  sim = 1
  seed = 1
end
mkpath(resultsdir("seir_cases", "generated_quantities"))
mkpath(resultsdir("seir_cases", "posterior_predictive"))


## Load Data
# choose sim 
if sim == 1
  all_dat = CSV.read("data/sim_data/scenario1_fitted_cases_obsdata.csv", DataFrame)
  overdisp_priors = CSV.read(datadir("sim_data", string("overdisp_priors_sim", sim, ".csv")), DataFrame)
  const phi_sd = overdisp_priors[1, :sd] 
  const phi_mean = overdisp_priors[1, :mean]
  ## Define Priors
  include(projectdir("src/prior_constants_seir_cases.jl"))
end 

dat = subset(all_dat, :seed => ByRow(x -> x == seed))
data_cases = dat[:, :total_cases]
obstimes = dat[:, :new_week]
param_change_times = obstimes[1:(end - 1)]



## Control Parameters
n_forecast_times = 0

## Define ODE
include(projectdir("src/seir_ode_log.jl"))

## Load Model
if fixed == true  
  include(projectdir("src/bayes_seir_cases_fixed.jl")) 

else 
  include(projectdir("src/bayes_seir_cases.jl")) 
end 

  

my_model = bayes_seir_cases(
  data_cases, 
  obstimes, 
  param_change_times, 
  true, 
  prob,
  1e-9,
  1e-6)


missing_cases = repeat([missing], length(data_cases))
  
my_model_forecast_missing = bayes_seir_cases(
  missing_cases,
  obstimes,
  param_change_times,
  true,
  prob,
  1e-9,
  1e-6)

if priors_only
    prior_samples = load(resultsdir("seir_cases",string("prior_samples_scenario", sim, ".jld2")))["prior_samples"]
    
    indices_to_keep = .!isnothing.(generated_quantities(my_model, prior_samples));
    
    prior_samples_randn = ChainsCustomIndex(prior_samples, indices_to_keep);
    
    
    Random.seed!(1234)
    
    
    prior_predictive_randn = predict(my_model_forecast_missing, prior_samples_randn)
    CSV.write(resultsdir("seir_cases", string("prior_predictive_scenario", sim, ".csv")), DataFrame(prior_predictive_randn))
    
    Random.seed!(1234)
    prior_gq_randn = get_gq_chains(my_model, prior_samples_randn);
    CSV.write(resultsdir("seir_cases", string("prior_generated_quantities_scenario", sim, ".csv")), DataFrame(prior_gq_randn))
    
    
    exit()
end

posterior_samples = load(resultsdir("seir_cases", "posterior_samples", string("posterior_samples", "_scenario", sim, "_seed", seed, ".jld2")))["posterior_samples"][:, :, 1:4]

indices_to_keep = .!isnothing.(generated_quantities(my_model, posterior_samples));

posterior_samples_randn = ChainsCustomIndex(posterior_samples, indices_to_keep);

Random.seed!(1234)
predictive_randn = predict(my_model_forecast_missing, posterior_samples_randn)
CSV.write(resultsdir("seir_cases", "posterior_predictive", string("posterior_predictive", "_scenario", sim, "_seed", seed, ".csv")), DataFrame(predictive_randn))

Random.seed!(1234)
gq_randn = get_gq_chains(my_model, posterior_samples_randn);
CSV.write(resultsdir("seir_cases", "generated_quantities", string("generated_quantities", "_scenario", sim, "_seed", seed, ".csv")), DataFrame(gq_randn))

posterior_df = DataFrame(posterior_samples)
CSV.write(resultsdir("seir_cases", "generated_quantities", string("posterior_df", "_scenario", sim, "_seed", seed, ".csv")), DataFrame(posterior_samples))

