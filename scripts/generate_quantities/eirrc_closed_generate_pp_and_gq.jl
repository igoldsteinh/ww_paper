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
### base scenario
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



### LA data
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

# 10 replicates--EIRR(10)
if sim == 3
  all_dat = CSV.read("data/sim_data/scenario1_fitted_genecount_obsdata.csv", DataFrame)
  dat = subset(all_dat, :seed => ByRow(x -> x == seed))
## Define Priors
include(projectdir("src/prior_constants_eirr_closed_scenario1.jl"))

subset_dat = dat[:, [:new_time, 
                     :log_gene_copies1, 
                     :log_gene_copies2, 
                     :log_gene_copies3, 
                     :log_gene_copies4, 
                     :log_gene_copies5, 
                     :log_gene_copies6, 
                     :log_gene_copies7,
                     :log_gene_copies8,
                     :log_gene_copies9,
                     :log_gene_copies10]]
long_dat = DataFrames.stack(subset_dat, [:log_gene_copies1, 
                                         :log_gene_copies2, 
                                         :log_gene_copies3, 
                                         :log_gene_copies4, 
                                         :log_gene_copies5, 
                                         :log_gene_copies6, 
                                         :log_gene_copies7,
                                         :log_gene_copies8,
                                         :log_gene_copies9,
                                         :log_gene_copies10])
long_dat = filter(:value => value -> value > 0, long_dat)
data_log_copies = long_dat[:, :value]
end 

### Mean of 3 replicates--EIRR(3 mean)
if sim == 4
  all_dat = CSV.read("data/sim_data/scenario1_fitted_genecount_obsdata.csv", DataFrame)
  dat = subset(all_dat, :seed => ByRow(x -> x == seed))
## Define Priors
include(projectdir("src/prior_constants_eirr_closed_scenario1.jl"))

subset_dat = dat[:, [:new_time, :log_mean_copiesthree]]
long_dat = subset_dat
data_log_copies = long_dat[:, :log_mean_copiesthree]
end 

### Mean of ten replicates--EIRR(10 mean)
if sim == 5
  all_dat = CSV.read("data/sim_data/scenario1_fitted_genecount_obsdata.csv", DataFrame)
  dat = subset(all_dat, :seed => ByRow(x -> x == seed))
## Define Priors
include(projectdir("src/prior_constants_eirr_closed_scenario1.jl"))

subset_dat = dat[:, [:new_time, :log_mean_copiesten]]
long_dat = subset_dat
data_log_copies = long_dat[:, :log_mean_copiesten]
end 

### 1 replicate--EIRR(1)
if sim == 6
  all_dat = CSV.read("data/sim_data/scenario1_fitted_genecount_obsdata.csv", DataFrame)
  dat = subset(all_dat, :seed => ByRow(x -> x == seed))
## Define Priors
include(projectdir("src/prior_constants_eirr_closed_scenario1.jl"))

subset_dat = dat[:, [:new_time, :log_gene_copies1]]
long_dat = subset_dat
data_log_copies = long_dat[:, :log_gene_copies1]
end 

### lambda centered at 0.8--Low Prop
if sim == 8
  all_dat = CSV.read("data/sim_data/scenario1_fitted_genecount_obsdata.csv", DataFrame)
  dat = subset(all_dat, :seed => ByRow(x -> x == seed))
## Define Priors
include(projectdir("src/prior_constants_eirr_closed_scenario8.jl"))

subset_dat = dat[:, [:new_time, :log_gene_copies1, :log_gene_copies2, :log_gene_copies3]]
long_dat = DataFrames.stack(subset_dat, [:log_gene_copies1, :log_gene_copies2, :log_gene_copies3])
long_dat = filter(:value => value -> value > 0, long_dat)
data_log_copies = long_dat[:, :value]

end 

### E and I initial centered low--Low Init
if sim == 9
  all_dat = CSV.read("data/sim_data/scenario1_fitted_genecount_obsdata.csv", DataFrame)
  dat = subset(all_dat, :seed => ByRow(x -> x == seed))
## Define Priors
include(projectdir("src/prior_constants_eirr_closed_scenario9.jl"))

subset_dat = dat[:, [:new_time, :log_gene_copies1, :log_gene_copies2, :log_gene_copies3]]
long_dat = DataFrames.stack(subset_dat, [:log_gene_copies1, :log_gene_copies2, :log_gene_copies3])
long_dat = filter(:value => value -> value > 0, long_dat)
data_log_copies = long_dat[:, :value]

end 

### E and I centered high--High Init
if sim == 10
  all_dat = CSV.read("data/sim_data/scenario1_fitted_genecount_obsdata.csv", DataFrame)
  dat = subset(all_dat, :seed => ByRow(x -> x == seed))
## Define Priors
include(projectdir("src/prior_constants_eirr_closed_scenario10.jl"))

subset_dat = dat[:, [:new_time, :log_gene_copies1, :log_gene_copies2, :log_gene_copies3]]
long_dat = DataFrames.stack(subset_dat, [:log_gene_copies1, :log_gene_copies2, :log_gene_copies3])
long_dat = filter(:value => value -> value > 0, long_dat)
data_log_copies = long_dat[:, :value]

end 

obstimes = long_dat[:, :new_time]
obstimes = convert(Vector{Float64}, obstimes)

# choose change times 
if maximum(obstimes) % 7 == 0
  param_change_max = maximum(obstimes) - 7
else 
  param_change_max = maximum(obstimes)
end 
param_change_times = collect(7:7.0:param_change_max)
full_time_series = collect(minimum(obstimes):grid_size:maximum(obstimes))
outs_tmp = dualcache(zeros(6,length(full_time_series)), 10)

index = zeros(length(obstimes))
for i in 1:length(index)
    time = obstimes[i]
    index[i] = indexin(time, full_time_series)[1]
end 

## Define closed form solution
include(projectdir("src/newnew_closed_soln_eirr_withincid.jl"))


## Load Model
include(projectdir("src/new_bayes_eirrc_closed.jl"))



my_model = new_bayes_eirrc_closed!(
    outs_tmp, 
    data_log_copies,
    obstimes, 
    param_change_times,
    grid_size,
    index)
  
  
  missing_log_copies = repeat([missing], length(data_log_copies))
  
  my_model_forecast_missing = new_bayes_eirrc_closed!(
    outs_tmp, 
    missing_log_copies,
    obstimes,
    param_change_times,
    grid_size,
    index)
  

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
