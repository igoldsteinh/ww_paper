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
using wastewater2
using PreallocationTools

sim =
if length(ARGS) == 0
   "real" 
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

# quick and dirty fix for doing this on hpc3
if sim == 56
  sim = "real"
end 
Logging.disable_logging(Logging.Warn)

mkpath(resultsdir("eirrc_closed"))
mkpath(resultsdir("eirrc_closed", "posterior_samples"))

## Control Parameters
n_samples = 250
n_chains = 4



## Load Data
# choose sim 
if sim == 1
  all_dat = CSV.read("data/sim_data/scenario1_fitted_genecount_obsdata.csv", DataFrame)
  dat = subset(all_dat, :seed => ByRow(x -> x == seed))
## Define Priors
include(projectdir("src/prior_constants_eirr_closed_scenario1.jl"))

subset_dat = dat[:, [:new_time, :log_gene_copies1, :log_gene_copies2, :log_gene_copies3]]
long_dat = DataFrames.stack(subset_dat, [:log_gene_copies1, :log_gene_copies2, :log_gene_copies3])
long_dat = filter(:value => value -> value > 0, long_dat)
# long_dat = filter(:value => value -> value < 16, long_dat)
data_log_copies = long_dat[:, :value]
end 




if sim == "real"
  dat = CSV.read("data/LA_daily_data_feb2022.csv", DataFrame)
  # for now I'm going to remove the last observation as we don't have a full week's worth of data for it
  dat = filter(:year_day => year_day -> year_day < 423, dat)
  ## Define Priors
include(projectdir("src", string("prior_constants_eirr_closed_LAdata", "_seed", seed, ".jl"))) 

subset_dat = dat[:, [:new_time, :log_gene_copies1, :log_gene_copies2, :log_gene_copies3]]
long_dat = DataFrames.stack(subset_dat, [:log_gene_copies1, :log_gene_copies2, :log_gene_copies3])
long_dat = filter(:value => value -> value > 0, long_dat)
data_log_copies = long_dat[:, :value]

end 

obstimes = long_dat[:, :new_time]
obstimes = convert(Vector{Float64}, obstimes)
# trying to avoid the stupid situation where we're telling to change at the end of the solver which doesn't make sense
if maximum(obstimes) % 7 == 0
  param_change_max = maximum(obstimes) - 7
else 
  param_change_max = maximum(obstimes)
end 
param_change_times = collect(7:7.0:param_change_max)
outs_tmp = dualcache(zeros(6,length(1:obstimes[end])), 10)

## Define closed form solution
include(projectdir("src/closed_soln_eirr_withincid.jl"))


## Load Model
include(projectdir("src/bayes_eirrc_closed.jl"))


my_model = bayes_eirrc_closed!(
    outs_tmp, 
    data_log_copies,
    obstimes, 
    param_change_times)
  
  
  missing_log_copies = repeat([missing], length(data_log_copies))
  
  my_model_forecast_missing = bayes_eirrc_closed!(
    outs_tmp, 
    missing_log_copies,
    obstimes,
    param_change_times)
  

if priors_only
    prior_samples = load(resultsdir("eirrc_closed",string("prior_samples_scenario", sim, "_seed", seed, ".jld2")))["prior_samples"]
    
    indices_to_keep = .!isnothing.(generated_quantities(my_model, prior_samples));
    
    prior_samples_randn = ChainsCustomIndex(prior_samples, indices_to_keep);
    
    
    Random.seed!(seed)
    
    
    prior_predictive_randn = predict(my_model_forecast_missing, prior_samples_randn)
    CSV.write(resultsdir("eirrc_closed", string("prior_predictive_scenario", sim, "_seed", seed, ".csv")), DataFrame(prior_predictive_randn))
    
    Random.seed!(seed)
    prior_gq_randn = get_gq_chains(my_model, prior_samples_randn);
    CSV.write(resultsdir("eirrc_closed", string("prior_generated_quantities_scenario", sim, "_seed", seed, ".csv")), DataFrame(prior_gq_randn))
    
    
        exit()
end

posterior_samples = load(resultsdir("eirrc_closed", "posterior_samples", string("posterior_samples", "_scenario", sim, "_seed", seed, ".jld2")))["posterior_samples"]

indices_to_keep = .!isnothing.(generated_quantities(my_model, posterior_samples));

posterior_samples_randn = ChainsCustomIndex(posterior_samples, indices_to_keep);


Random.seed!(seed)
predictive_randn = predict(my_model_forecast_missing, posterior_samples_randn)
CSV.write(resultsdir("eirrc_closed", "posterior_predictive", string("posterior_predictive", "_scenario", sim, "_seed", seed,  ".csv")), DataFrame(predictive_randn))

Random.seed!(seed)
gq_randn = get_gq_chains(my_model, posterior_samples_randn);
CSV.write(resultsdir("eirrc_closed", "generated_quantities", string("generated_quantities", "_scenario", sim, "_seed", seed, ".csv")), DataFrame(gq_randn))


posterior_df = DataFrame(posterior_samples)
CSV.write(resultsdir("eirrc_closed", "generated_quantities", string("posterior_df", "_scenario", sim, "_seed", seed, ".csv")), DataFrame(posterior_samples))
