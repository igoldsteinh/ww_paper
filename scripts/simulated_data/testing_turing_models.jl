# testing models with different closed solutions 
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
using BenchmarkTools

#just use first data set 
sim = 1
seed = 1
all_dat = CSV.read("data/sim_data/scenario1_fitted_genecount_obsdata.csv", DataFrame)
dat = subset(all_dat, :seed => ByRow(x -> x == seed))
## Define Priors
include(projectdir("src/prior_constants_eirr_closed_scenario1.jl"))

subset_dat = dat[:, [:new_time, :log_gene_copies1, :log_gene_copies2, :log_gene_copies3]]
long_dat = DataFrames.stack(subset_dat, [:log_gene_copies1, :log_gene_copies2, :log_gene_copies3])
long_dat = filter(:value => value -> value > 0, long_dat)
# long_dat = filter(:value => value -> value < 16, long_dat)
data_log_copies = long_dat[:, :value]
grid_size = 1.0

n_samples = 20
n_chains = 1


obstimes = long_dat[:, :new_time]
obstimes = convert(Vector{Float64}, obstimes)

# pick the change times 
if maximum(obstimes) % 7 == 0
  param_change_max = maximum(obstimes) - 7
else 
  param_change_max = maximum(obstimes)
end 
param_change_times = collect(7:7.0:param_change_max)
full_time_series = collect(unique(obstimes)[1]:grid_size:unique(obstimes)[end])
outs_tmp = dualcache(zeros(6,length(full_time_series)), 10)

index = zeros(length(obstimes))
for i in 1:length(index)
    time = obstimes[i]
    index[i] = indexin(time, full_time_series)[1]
end 
# make old model 
## Define closed form solution
include(projectdir("src/closed_soln_eirr_withincid.jl"))


## Load Model
include(projectdir("src/bayes_eirrc_closed.jl"))


my_model = bayes_eirrc_closed!(
    outs_tmp, 
    data_log_copies,
    obstimes, 
    param_change_times)


Random.seed!(seed)

MAP_init = optimize_many_MAP(my_model, 10, 1, true)[1]
    
Random.seed!(seed)
MAP_noise = vcat(randn(length(MAP_init) - 1, n_chains), transpose(zeros(n_chains)))
MAP_noise = [MAP_noise[:,i] for i in 1:size(MAP_noise,2)]
    
init = repeat([MAP_init], n_chains) .+ 0.05 * MAP_noise
    
# make new model 
## Define closed form solution
include(projectdir("src/newnew_closed_soln_eirr_withincid.jl"))


## Load Model
include(projectdir("src/new_bayes_eirrc_closed.jl"))


my_model_new = new_bayes_eirrc_closed!(
    outs_tmp, 
    data_log_copies,
    obstimes, 
    param_change_times,
    grid_size, 
    index)


Random.seed!(seed)

MAP_init = optimize_many_MAP(my_model_new, 10, 1, true)[1]
    
Random.seed!(seed)
MAP_noise = vcat(randn(length(MAP_init) - 1, n_chains), transpose(zeros(n_chains)))
MAP_noise = [MAP_noise[:,i] for i in 1:size(MAP_noise,2)]
    
init_new = repeat([MAP_init], n_chains) .+ 0.05 * MAP_noise

# fit each model 
Random.seed!(seed)
@btime posterior_samples = sample(my_model, NUTS(), MCMCThreads(), n_samples, n_chains, discard_initial = n_samples, init_params = init)
@btime posterior_samples_new = sample(my_model_new, NUTS(), MCMCThreads(), n_samples, n_chains, discard_initial = n_samples, init_params = init_new)