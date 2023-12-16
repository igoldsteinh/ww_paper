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
using LineSearches
using ww_paper
using Logging

sim =
if length(ARGS) == 0
  1
else
  parse(Int64, ARGS[1])
end

priors_only = sim == 0

seed = 
if length(ARGS) == 0
  2
else 
  parse(Int64, ARGS[2])
end 


if priors_only
  sim = 1
  seed = 1
end


  
  
Logging.disable_logging(Logging.Warn)
mkpath(resultsdir("seir_cases"))
mkpath(resultsdir("seir_cases", "posterior_samples"))

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
obstimes = convert(Vector{Float64}, obstimes)
param_change_times = obstimes[1:(end - 1)]

## Control Parameters
n_samples = 500
n_chains = 4


## Define ODE
include(projectdir("src/seir_ode_log.jl"))

## Load Model
if fixed == true  
  include(projectdir("src/bayes_seir_cases_fixed.jl")) 

else 
  include(projectdir("src/bayes_seir_cases.jl")) 
end 

  
    

my_model_optimize = bayes_seir_cases(
    data_cases, 
    obstimes, 
    param_change_times, 
    true, 
    prob,
    1e-11,
    1e-8)

my_model = bayes_seir_cases(
  data_cases, 
  obstimes, 
  param_change_times, 
  true, 
  prob,
  1e-9,
  1e-6)


# Sample Posterior

if priors_only
  Random.seed!(seed)
  prior_samples = sample(my_model, Prior(), MCMCThreads(), 400, n_chains)
  wsave(resultsdir("seir_cases", string("prior_samples_scenario", sim, ".jld2")), @dict prior_samples)
  exit()
end
Random.seed!(seed)
MAP_init = optimize_many_MAP(my_model_optimize, 10, 1, true)[1]

Random.seed!(seed)
MAP_noise = randn(length(MAP_init), n_chains)
MAP_noise = [MAP_noise[:,i] for i in 1:size(MAP_noise,2)]

init = repeat([MAP_init], n_chains) .+ 0.05 * MAP_noise
Random.seed!(seed)

posterior_samples = sample(my_model, NUTS(-1, 0.8), MCMCThreads(), n_samples, n_chains, discard_initial = n_samples, init_params = init)



wsave(resultsdir("seir_cases", "posterior_samples", string("posterior_samples_scenario", sim, "_seed", seed, ".jld2")), @dict posterior_samples)
