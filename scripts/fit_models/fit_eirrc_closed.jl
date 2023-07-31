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
using wastewater2
using Logging
using PreallocationTools

sim =
if length(ARGS) == 0
   "ODE"
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

### Baseline scenario
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


### Baseline scenario
if sim == "ODE"
  all_dat = CSV.read("data/sim_data/ODE_comp_fitted_genecount_obsdata.csv", DataFrame)
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

### 10 replicates--EIRR (10)
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

### mean of 3 replicates--EIRR (3 mean)
if sim == 4
  all_dat = CSV.read("data/sim_data/scenario1_fitted_genecount_obsdata.csv", DataFrame)
  dat = subset(all_dat, :seed => ByRow(x -> x == seed))
## Define Priors
include(projectdir("src/prior_constants_eirr_closed_scenario1.jl"))

subset_dat = dat[:, [:new_time, :log_mean_copiesthree]]
long_dat = subset_dat
data_log_copies = long_dat[:, :log_mean_copiesthree]
end 

### mean of 10 replicates-- EIRR (10 mean)
if sim == 5
  all_dat = CSV.read("data/sim_data/scenario1_fitted_genecount_obsdata.csv", DataFrame)
  dat = subset(all_dat, :seed => ByRow(x -> x == seed))
## Define Priors
include(projectdir("src/prior_constants_eirr_closed_scenario1.jl"))

subset_dat = dat[:, [:new_time, :log_mean_copiesten]]
long_dat = subset_dat
data_log_copies = long_dat[:, :log_mean_copiesten]
end 

### 1 replicate--EIRR (1)
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

### E and I iniital compartments centered low--Low Init
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

### E and I iniital compartments centered high--High Init
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

# pick the change times 
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

# Sample Posterior

if priors_only
  Random.seed!(seed)
  prior_samples = sample(my_model, Prior(), MCMCThreads(), 400, n_chains)
  wsave(resultsdir("eirrc_closed",string("prior_samples_scenario", sim, "_seed", seed, ".jld2")), @dict prior_samples)
  exit()
end

Random.seed!(seed)

MAP_init = optimize_many_MAP(my_model, 10, 1, true)[1]

Random.seed!(seed)
MAP_noise = vcat(randn(length(MAP_init) - 1, n_chains), transpose(zeros(n_chains)))
MAP_noise = [MAP_noise[:,i] for i in 1:size(MAP_noise,2)]

init = repeat([MAP_init], n_chains) .+ 0.05 * MAP_noise

Random.seed!(seed)
posterior_samples = sample(my_model, NUTS(), MCMCThreads(), n_samples, 4, discard_initial = n_samples, init_params = init)
wsave(resultsdir("eirrc_closed", "posterior_samples", string("posterior_samples_scenario", sim, "_seed", seed, ".jld2")), @dict posterior_samples)


