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
using PreallocationTools


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
mkpath(resultsdir("eir_cases", "generated_quantities"))
mkpath(resultsdir("eir_cases", "posterior_predictive"))

## Load Data
# choose sim 
if sim == 1
  all_dat = CSV.read("data/sim_data/scenario1_fitted_cases_obsdata.csv", DataFrame)
  overdisp_priors = CSV.read(datadir("sim_data", string("overdisp_priors_sim", sim, ".csv")), DataFrame)
  const phi_sd = overdisp_priors[1, :sd] 
  const phi_mean = overdisp_priors[1, :mean]
  dat = subset(all_dat, :seed => ByRow(x -> x == seed))

  ## Define Priors
  include(projectdir("src/prior_constants_eir_cases.jl"))
end 

if sim == "real"
  all_dat = CSV.read("data/LA_EIR_data.csv", DataFrame)
  dat = all_dat
  overdisp_priors = CSV.read(datadir("sim_data", string("overdisp_priors_sim", sim, ".csv")), DataFrame)
  const phi_sd = overdisp_priors[1, :sd] 
  const phi_mean = overdisp_priors[1, :mean]

  ## Define Priors
  include(projectdir("src/prior_constants_eir_cases_LA.jl"))
end 


data_cases = dat[:, :total_cases]
obstimes = dat[:, :new_week]
obstimes = convert(Vector{Float64}, obstimes)
param_change_times = obstimes[1:(end - 1)]

outs_tmp = dualcache(zeros(5,length(1:obstimes[end])), 10)


## Control Parameters
n_forecast_times = 0

## Define Priors
include(projectdir("src/prior_constants_eir_cases.jl"))

## Define closed form solution
include(projectdir("src/closed_soln_eir.jl"))


## Load Model
include(projectdir("src/bayes_eir_cases.jl"))

my_model = bayes_eir_cases!(
    outs_tmp, 
    data_cases, 
    obstimes, 
    param_change_times)


  missing_cases = repeat([missing], length(data_cases))
  
  my_model_forecast_missing = bayes_eir_cases!(
    outs_tmp, 
    missing_cases, 
    obstimes, 
    param_change_times)

if priors_only
    prior_samples = load(resultsdir("eir_cases",string("prior_samples_scenario", sim, "_seed", seed, ".jld2")))["prior_samples"]
    
    indices_to_keep = .!isnothing.(generated_quantities(my_model, prior_samples));
    
    prior_samples_randn = ChainsCustomIndex(prior_samples, indices_to_keep);
    
    
    Random.seed!(seed)
    
    
    
    Random.seed!(seed)
    prior_gq_randn = get_gq_chains(my_model, prior_samples_randn);
    CSV.write(resultsdir("eir_cases", string("prior_generated_quantities_scenario", sim, "_seed", seed, ".csv")), DataFrame(prior_gq_randn))
    
    
    exit()
end

posterior_samples = load(resultsdir("eir_cases", "posterior_samples", string("posterior_samples_scenario", sim, "_seed", seed, ".jld2")))["posterior_samples"][:, :, 1:4]

indices_to_keep = .!isnothing.(generated_quantities(my_model, posterior_samples));

posterior_samples_randn = ChainsCustomIndex(posterior_samples, indices_to_keep);



Random.seed!(seed)
gq_randn = get_gq_chains(my_model, posterior_samples_randn);
CSV.write(resultsdir("eir_cases", "generated_quantities", string("generated_quantities", "_scenario", sim, "_seed", seed, ".csv")), DataFrame(gq_randn))

posterior_df = DataFrame(posterior_samples)
CSV.write(resultsdir("eir_cases", "generated_quantities", string("posterior_df", "_scenario", sim, "_seed", seed,  ".csv")), DataFrame(posterior_samples))

# for simulations I am going to not record posterior predictives
# this can be done easily if needed by uncommenting the code below
# the reason to not automate this is that turing begins by plugging in arbitrary values into the model before using the values in the actual chain 
# these arbitrary values can cause numerical errors which cause the whole script to fail 
# to fix the problem, just change the seed, this will change the initial values turing uses 

# Random.seed!(seed)
# predictive_randn = predict(my_model_forecast_missing, posterior_samples_randn)
# CSV.write(resultsdir("eir_cases", "posterior_predictive", string("posterior_predictive", "_scenario", sim, "_seed", seed, ".csv")), DataFrame(predictive_randn))
