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

seed = 
if length(ARGS) == 0
  1
else 
  parse(Int64, ARGS[2])
end 



priors_only = sim == 0

if priors_only
  sim = 1
  seed = 1
end

mkpath(resultsdir("seirr_student", "generated_quantities"))
mkpath(resultsdir("seirr_student", "posterior_predictive"))


## Load Data
# choose sim 
if sim == 1
  all_dat = CSV.read("data/sim_data/scenario1_fitted_genecount_obsdata.csv", DataFrame)
  dat = subset(all_dat, :seed => ByRow(x -> x == seed))
  subset_dat = dat[:, [:new_time, :log_gene_copies1, :log_gene_copies2, :log_gene_copies3]]
  long_dat = DataFrames.stack(subset_dat, [:log_gene_copies1, :log_gene_copies2, :log_gene_copies3])
  data_log_copies = long_dat[:, :value]
end 

obstimes = long_dat[:, :new_time]
obstimes = convert(Vector{Float64}, obstimes)
if maximum(obstimes) % 7 == 0
  param_change_max = maximum(obstimes) - 7
else 
  param_change_max = maximum(obstimes)
end 
param_change_times = collect(7:7.0:param_change_max)
  
## Define Priors
include(projectdir("src/prior_constants_seirr_student.jl"))

## Define ODE
include(projectdir("src/seirr_ode_log.jl"))

## Load Model
include(projectdir("src/bayes_seirr_student.jl"))
  
my_model = bayes_seirr_student(
  data_log_copies,
  obstimes, 
  param_change_times, 
  true,
  prob,
  1e-9,
  1e-6)

missing_log_copies = repeat([missing], length(data_log_copies))
  
my_model_forecast_missing = bayes_seirr_student(
  missing_log_copies,
  obstimes,
  param_change_times,
  true,
  prob,
  1e-9,
  1e-6)
  

if priors_only
    prior_samples = load(resultsdir("seirr_student", string("prior_samples_scenario", sim, ".jld2")))["prior_samples"]
    
    indices_to_keep = .!isnothing.(generated_quantities(my_model, prior_samples));
    
    prior_samples_randn = ChainsCustomIndex(prior_samples, indices_to_keep);
    
    
    Random.seed!(seed)
    
    
    prior_predictive_randn = predict(my_model_forecast_missing, prior_samples_randn)
    CSV.write(resultsdir("seirr_student", string("prior_predictive_scenario", sim, ".csv")), DataFrame(prior_predictive_randn))
    
    Random.seed!(seed)
    prior_gq_randn = get_gq_chains(my_model, prior_samples_randn);
    CSV.write(resultsdir("seirr_student", string("prior_generated_quantities_scenario", sim, ".csv")), DataFrame(prior_gq_randn))
    
    
    exit()
end

posterior_samples = load(resultsdir("seirr_student", "posterior_samples", string("posterior_samples", "_scenario", sim, "_seed", seed, ".jld2")))["posterior_samples"]

indices_to_keep = .!isnothing.(generated_quantities(my_model, posterior_samples));

posterior_samples_randn = ChainsCustomIndex(posterior_samples, indices_to_keep);

Random.seed!(seed)
predictive_randn = predict(my_model_forecast_missing, posterior_samples_randn)
CSV.write(resultsdir("seirr_student", "posterior_predictive", string("posterior_predictive", "_scenario", sim, "_seed", seed,  ".csv")), DataFrame(predictive_randn))

Random.seed!(seed)
gq_randn = get_gq_chains(my_model, posterior_samples_randn);
CSV.write(resultsdir("seirr_student", "generated_quantities", string("generated_quantities", "_scenario", sim, "_seed", seed, ".csv")), DataFrame(gq_randn))

posterior_df = DataFrame(posterior_samples)
CSV.write(resultsdir("seirr_student", "generated_quantities", string("posterior_df", "_scenario", sim, "_seed", seed, ".csv")), DataFrame(posterior_samples))

