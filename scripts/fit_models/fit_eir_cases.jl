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
using PreallocationTools

sim =
if length(ARGS) == 0
  "real"
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

# quick and dirty fix for doing this on hpc3
if sim == 56
  sim = "real"
end 


Logging.disable_logging(Logging.Warn)
mkpath(resultsdir("eir_cases"))
mkpath(resultsdir("eir_cases", "posterior_samples"))

## Control Parameters
n_samples = 500
n_chains = 4


## Define closed form solution
include(projectdir("src/closed_soln_eir.jl"))

## Load Model
include(projectdir("src/bayes_eir_cases.jl"))

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
param_change_times = convert(Vector{Float64}, param_change_times)

outs_tmp = dualcache(zeros(5,length(1:obstimes[end])), 10)

    

my_model = bayes_eir_cases!(
    outs_tmp, 
    data_cases, 
    obstimes, 
    param_change_times)


if priors_only
  Random.seed!(seed)
  prior_samples = sample(my_model, Prior(), MCMCThreads(), 400, n_chains)
  wsave(resultsdir("eir_cases", string("prior_samples_scenario", sim, "_seed", seed, ".jld2")), @dict prior_samples)
  exit()
end
Random.seed!(seed)
MAP_init = optimize_many_MAP(my_model, 10, 1, true)[1]

Random.seed!(seed)
MAP_noise = randn(length(MAP_init), n_chains)
MAP_noise = [MAP_noise[:,i] for i in 1:size(MAP_noise,2)]

init = repeat([MAP_init], n_chains) .+ 0.05 * MAP_noise
Random.seed!(seed)

posterior_samples = sample(my_model, NUTS(-1, 0.8), MCMCThreads(), n_samples, n_chains, discard_initial = n_samples, init_params = init)


wsave(resultsdir("eir_cases", "posterior_samples", string("posterior_samples_scenario", sim, "_seed", seed, ".jld2")), @dict posterior_samples)
