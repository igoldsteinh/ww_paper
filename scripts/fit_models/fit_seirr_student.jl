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
102
else
  parse(Int64, ARGS[1])
end

seed = 
if length(ARGS) == 0
  1
else 
  parse(Int64, ARGS[2])
end 

if sim == "fixed" || sim == "apf"
fixed = true 
else 
fixed = false 
end 
priors_only = sim == 0

if sim == "frw"
  fixed_rw = true
else
  fixed_rw = false 
end 

if priors_only
  sim = 1
  seed = 1
end

Logging.disable_logging(Logging.Warn)
mkpath(resultsdir("seirr_student"))
mkpath(resultsdir("seirr_student", "posterior_samples"))

## Control Parameters
n_samples = 500
n_chains = 4

## Load Data
# choose sim 
if sim == 1
  all_dat = CSV.read("data/sim_data/scenario1_fitted_genecount_obsdata.csv", DataFrame)
  dat = subset(all_dat, :seed => ByRow(x -> x == seed))
  subset_dat = dat[:, [:new_time, :log_gene_copies1, :log_gene_copies2, :log_gene_copies3]]
  long_dat = DataFrames.stack(subset_dat, [:log_gene_copies1, :log_gene_copies2, :log_gene_copies3])
  data_log_copies = long_dat[:, :value]

  ## Define Priors
  include(projectdir("src/prior_constants_seirr_student.jl"))
end 

if sim == 102
  all_dat = CSV.read("data/sim_data/scenario1_fitted_genecount_obsdata.csv", DataFrame)
  dat = subset(all_dat, :seed => ByRow(x -> x == seed))
  subset_dat = dat[:, [:new_time, :log_gene_copies1, :log_gene_copies2, :log_gene_copies3]]
  long_dat = DataFrames.stack(subset_dat, [:log_gene_copies1, :log_gene_copies2, :log_gene_copies3])
  data_log_copies = long_dat[:, :value]

  ## Define Priors
  include(projectdir("src/prior_constants_seirr_student_scenario102.jl"))
end 



obstimes = long_dat[:, :new_time]
obstimes = convert(Vector{Float64}, obstimes)
if sim != 102
  if maximum(obstimes) % 7 == 0
    param_change_max = maximum(obstimes) - 7
  else 
    param_change_max = maximum(obstimes)
  end 
  param_change_times = collect(7:7.0:param_change_max)
else 
  if maximum(obstimes) % 3 == 0
    param_change_max = maximum(obstimes) - 3
  else 
    param_change_max = maximum(obstimes)
  end 
  param_change_times = collect(3:3.0:param_change_max)
end 

  
  

## Define ODE
include(projectdir("src/seirr_ode_log.jl"))

## Load Model
if fixed == true
include(projectdir("src/bayes_seirr_student_fixed.jl"))
elseif fixed_rw == true 
include(projectdir("src/bayes_seirr_student_frw.jl"))
else 
include(projectdir("src/bayes_seirr_student.jl"))
end

if sim == 102
  my_model_optimize = bayes_seirr_student(
      data_log_copies, 
      obstimes, 
      param_change_times, 
      true, 
      prob,
      1e-11,
      1e-8,
      3)

  my_model = bayes_seirr_student(
    data_log_copies, 
    obstimes, 
    param_change_times, 
    true, 
    prob,
    1e-9,
    1e-6, 
    3)
else 
  my_model_optimize = bayes_seirr_student(
      data_log_copies, 
      obstimes, 
      param_change_times, 
      true, 
      prob,
      1e-11,
      1e-8,
      7)

  my_model = bayes_seirr_student(
    data_log_copies, 
    obstimes, 
    param_change_times, 
    true, 
    prob,
    1e-9,
    1e-6,
    7)
end


# Sample Posterior

if priors_only
  Random.seed!(seed)
  prior_samples = sample(my_model, Prior(), MCMCThreads(), n_samples, n_chains)
  wsave(resultsdir("seirr_student", string("prior_samples_scenario", sim, ".jld2")), @dict prior_samples)
  exit()
end

Random.seed!(seed)
MAP_init = optimize_many_MAP(my_model, 10, 1, true)[1]

Random.seed!(seed)
# remove noise for the last parameter df, which is always positive and should not become negative due to noise
# sub it with zeros instead 
MAP_noise = vcat(randn(length(MAP_init) - 1, n_chains), transpose(zeros(n_chains)))
MAP_noise = [MAP_noise[:,i] for i in 1:size(MAP_noise,2)]

init = repeat([MAP_init], n_chains) .+ 0.05 * MAP_noise

Random.seed!(seed)

posterior_samples = sample(my_model, NUTS(), MCMCThreads(), n_samples, n_chains, discard_initial = n_samples, init_params = init)

wsave(resultsdir("seirr_student", "posterior_samples", string("posterior_samples_scenario", sim, "_seed", seed, ".jld2")), @dict posterior_samples)
