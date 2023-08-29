# testing out our new function which allows for finer grids
using DifferentialEquations
using DrWatson
using BenchmarkTools
using CSV
using DataFrames
include(projectdir("src/eirr_ode_log.jl"))
include(projectdir("src/closed_soln_eirr_withincid.jl"))
include(projectdir("src/newnew_closed_soln_eirr_withincid.jl"))
include(projectdir("src/eirr_ode.jl"))

E_init = 5
I_init = 5
R1_init = 5
R2_init = 5
C_init = 1
u0 = [E_init, I_init, R1_init, R2_init, C_init]
gamma = 1/4
nu = 1/7
eta = 1/18
rt = [1.5, 1.25]
alphas = rt .* nu 
times = collect(1:1:14)
param_change_times = 7.0
t0 = 0
grid_size = 1
outs_tmp = dualcache(zeros(6,length(unique(times))), 10)
init_conds = u0
change_times = param_change_times

my_sol = new_eirrc_closed_solution!(outs_tmp, times, param_change_times, t0, alphas, init_conds, gamma, nu, eta)

new_my_sol = newnew_eirrc_closed_solution!(outs_tmp, times, param_change_times, grid_size, t0, alphas, init_conds, gamma, nu, eta)


times =[1:1.0:14.0;]
change_times = [7.0]
t0 = 0.0
alphas_no_init = alphas[2]

log_prob = ODEProblem(eirr_ode_log!,
log.(u0),
(0.0, times[end]),
[alphas[1], gamma, nu, eta])

prob = ODEProblem(eirr_ode!,
u0,
(0.0, times[end]),
[alphas[1], gamma, nu, eta])

function param_affect_β_IFR!(integrator)
ind_t = searchsortedfirst(change_times, integrator.t) # Find the index of param_change_times that contains the current timestep
integrator.p[1] = alphas_no_init[ind_t] # Replace β with a new value from β_t_values
end

param_callback = PresetTimeCallback(change_times, param_affect_β_IFR!, save_positions = (false, false))

# Solve the ODE at intervals of 1.0, could also solve at obstimes
log_sol = solve(log_prob, Tsit5(), callback = param_callback, saveat = 1.0, save_start = true, verbose = false, abstol = 1e-9, reltol = 1e-6)  
ode_sol = solve(prob, Tsit5(), callback = param_callback, saveat = 1.0, save_start = true, verbose = false, abstol = 1e-9, reltol = 1e-6)  
exp_sol = exp.(sol)
# eir_closed_solution!(outs_tmp, times, change_times, 0.0, alphas, init_conds, gamma, nu)

### testing on time series length of scenario 1
all_dat = CSV.read("data/sim_data/scenario1_fitted_genecount_obsdata.csv", DataFrame)
seed = 1
dat = subset(all_dat, :seed => ByRow(x -> x == seed))

subset_dat = dat[:, [:new_time, :log_gene_copies1, :log_gene_copies2, :log_gene_copies3]]
long_dat = DataFrames.stack(subset_dat, [:log_gene_copies1, :log_gene_copies2, :log_gene_copies3])
long_dat = filter(:value => value -> value > 0, long_dat)
# long_dat = filter(:value => value -> value < 16, long_dat)
data_log_copies = long_dat[:, :value]
grid_size = 1.0


obstimes = long_dat[:, :new_time]
obstimes = convert(Vector{Float64}, obstimes)

# pick the change times 
if maximum(obstimes) % 7 == 0
  param_change_max = maximum(obstimes) - 7
else 
  param_change_max = maximum(obstimes)
end 
param_change_times = collect(7:7.0:param_change_max)
param_change_times = collect(7:7.0:param_change_max)
full_time_series = collect(minimum(obstimes):grid_size:maximum(obstimes))
outs_tmp = dualcache(zeros(6,length(full_time_series)), 10)

E_init = 5
I_init = 5
R1_init = 5
R2_init = 5
C_init = 1
u0 = [E_init, I_init, R1_init, R2_init, C_init]
gamma = 1/4
nu = 1/7
eta = 1/18
rt = collect(range(1.1, 1.5, length(param_change_times) + 1))
alphas = rt .* nu 
times = full_time_series 
t0 = 0
change_times = param_change_times
init_conds = u0
@btime my_sol = newnew_eirrc_closed_solution!(outs_tmp, times, param_change_times, grid_size, t0, alphas, init_conds, gamma, nu, eta)


alphas_no_init = alphas[2:end]

log_prob = ODEProblem(eirr_ode_log!,
log.(u0),
(0.0, maximum(times)),
[alphas[1], gamma, nu, eta])

prob = ODEProblem(eirr_ode!,
u0,
(0.0, maximum(times)),
[alphas[1], gamma, nu, eta])

function param_affect_β_IFR!(integrator)
ind_t = searchsortedfirst(change_times, integrator.t) # Find the index of param_change_times that contains the current timestep
integrator.p[1] = alphas_no_init[ind_t] # Replace β with a new value from β_t_values
end

param_callback = PresetTimeCallback(change_times, param_affect_β_IFR!, save_positions = (false, false))

# Solve the ODE at intervals of 1.0, could also solve at obstimes
@btime log_sol = solve(log_prob, Tsit5(), callback = param_callback, saveat = 1.0, save_start = true, verbose = false, abstol = 1e-11, reltol = 1e-8)  

# this did not help at all in fact it was even faster than the log version 
@btime ode_sol = solve(prob, Tsit5(), callback = param_callback, saveat = 1.0, save_start = true, verbose = false, abstol = 1e-9, reltol = 1e-6)  


# testing purely the matrix exponential 
alpha = alphas[1]
@btime my_fixed_alpha = fixed_alpha_closed_soln(outs_tmp, times,  grid_size, t0, alpha, init_conds, gamma, nu, eta)

prob = ODEProblem(eirr_ode!,
u0,
(0.0, maximum(times)),
[alphas[1], gamma, nu, eta])

@btime ode_fixed_alpha = solve(prob, Tsit5(), saveat = 1.0, save_start = true, verbose = false, abstol = 1e-9, reltol = 1e-6)  


log_prob = ODEProblem(eirr_ode_log!,
log.(u0),
(0.0, maximum(times)),
[alphas[1], gamma, nu, eta])

@btime log_ode_fixed_alpha = solve(log_prob, Tsit5(), saveat = 1.0, save_start = true, verbose = false, abstol = 1e-11, reltol = 1e-8)  

# testing speed on time series equal to length of real data
dat = CSV.read("data/LA_daily_data_feb2022.csv", DataFrame)
# for now I'm going to remove the last observation as we don't have a full week's worth of data for it
dat = filter(:year_day => year_day -> year_day < 423, dat)
## Define Priors
include(projectdir("src", string("prior_constants_eirr_closed_LAdata", "_seed", seed, ".jl"))) 

subset_dat = dat[:, [:new_time, :log_gene_copies1, :log_gene_copies2, :log_gene_copies3]]
long_dat = DataFrames.stack(subset_dat, [:log_gene_copies1, :log_gene_copies2, :log_gene_copies3])
long_dat = filter(:value => value -> value > 0, long_dat)
data_log_copies = long_dat[:, :value]
grid_size = 1.0
obstimes = long_dat[:, :new_time]
obstimes = convert(Vector{Float64}, obstimes)

# pick the change times 
if maximum(obstimes) % 7 == 0
  param_change_max = maximum(obstimes) - 7
else 
  param_change_max = maximum(obstimes)
end 
param_change_times = collect(7:7.0:param_change_max)
param_change_times = collect(7:7.0:param_change_max)
full_time_series = collect(minimum(obstimes):grid_size:maximum(obstimes))
outs_tmp = dualcache(zeros(6,length(full_time_series)), 10)

E_init = 5
I_init = 5
R1_init = 5
R2_init = 5
C_init = 1
u0 = [E_init, I_init, R1_init, R2_init, C_init]
gamma = 1/4
nu = 1/7
eta = 1/18
rt = fill(1.5, length(param_change_times) + 1)
alphas = rt .* nu 
times = obstimes 
t0 = 0
change_times = param_change_times
init_conds = u0
@btime my_sol = newnew_eirrc_closed_solution!(outs_tmp, times, param_change_times, grid_size, t0, alphas, init_conds, gamma, nu, eta)



alphas_no_init = alphas[2:end]
gamma = 1/4
nu = 1/7
eta = 1/18

prob = ODEProblem(eirr_ode_log!,
log.(u0),
(0.0, maximum(times)),
[alphas[1], gamma, nu, eta])

function param_affect_β_IFR!(integrator)
ind_t = searchsortedfirst(change_times, integrator.t) # Find the index of param_change_times that contains the current timestep
integrator.p[1] = alphas_no_init[ind_t] # Replace β with a new value from β_t_values
end

param_callback = PresetTimeCallback(change_times, param_affect_β_IFR!, save_positions = (false, false))

# Solve the ODE at intervals of 1.0, could also solve at obstimes
@btime sol = solve(prob, Tsit5(), callback = param_callback, saveat = 1.0, save_start = true, verbose = false, abstol = 1e-9, reltol = 1e-6)  
