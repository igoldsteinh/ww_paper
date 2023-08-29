# Closed form solution for EIRRC model with flexible resolution (daily, half-daily etc.)
using ForwardDiff
using PreallocationTools

function power(a,b)
  a^b
end 

# rewrite function for gradient friendly stuff
function correct_eirrc_matrix_exp(t, x::Array{T}) where {T<:Real}
  (alpha, gamma, nu, eta) = x
  matrix_exp = Array{Real}(undef, 5,5)
  # Row 1
matrix_exp[1,1] = -0.25*(exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
  (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
  (-2*alpha*gamma + gamma*nu - power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
(gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) - 
(exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
 (-gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
 (2*alpha*gamma - gamma*nu + power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
(4*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))
matrix_exp[1,2] = (alpha*exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
 (-gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
 (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
(4*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) + 
(alpha*exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
 (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
 (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
(4*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))
matrix_exp[1,3] = 0
matrix_exp[1,4] = 0
matrix_exp[1,5] = 0
# Row 2
matrix_exp[2,1] = -0.25*(exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
  (-2*alpha*gamma + gamma*nu - power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
  (2*alpha*gamma - gamma*nu + power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
(alpha*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) - 
(exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
 (2*alpha*gamma - gamma*nu + power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
 (2*alpha*gamma - gamma*nu + power(nu,2) + nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
(4*alpha*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))
matrix_exp[2,2] = (exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
 (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
 (2*alpha*gamma - gamma*nu + power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
(4*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) + 
(exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
 (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
 (2*alpha*gamma - gamma*nu + power(nu,2) + nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
(4*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))
matrix_exp[2,3] = 0
matrix_exp[2,4] = 0
matrix_exp[2,5] = 0
# Row 3
matrix_exp[3,1] = -0.25*(exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*nu*
  (gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
  (gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
  (-2*alpha*gamma + gamma*nu - power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
(alpha*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
  (-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
(exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*nu*
 (gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
 (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
 (2*alpha*gamma - gamma*nu + power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
(4*alpha*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
 (-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
((2*nu)/(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) - 
 (2*power(nu,2))/(alpha*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
 (2*nu)/(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) + 
 (2*power(nu,2))/(alpha*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))))/
(exp(eta*t)*(-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
   (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma)))

matrix_exp[3,2] = (exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*nu*
 (gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
 (gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
 (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
(4*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
 (-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
(exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*nu*
 (gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
 (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
 (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
(4*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
 (-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
((nu*(-gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
    (gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
    (gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
  (4*alpha*power(gamma,2)*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
 (nu*(gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
    (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
    (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
  (4*alpha*power(gamma,2)*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))))/
(exp(eta*t)*(-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
   (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma)))

matrix_exp[3,3] = -((sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma - 
  (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma))/
(exp(eta*t)*(-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
    (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma))))
matrix_exp[3,4] = 0
matrix_exp[3,5] = 0
# Row 4
matrix_exp[4,1] = -0.5*(exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*eta*nu*
  (-gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
  (-2*alpha*gamma + gamma*nu - power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
(alpha*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
  (-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
(exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*eta*nu*
 (-gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
 (2*alpha*gamma - gamma*nu + power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
(2*alpha*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
 (-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
((2*nu)/(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) - 
 (2*power(nu,2))/(alpha*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
 (2*nu)/(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) + 
 (2*power(nu,2))/(alpha*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))))/
(exp(eta*t)*(-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
   (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma))) + 
((-2*nu)/(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) - 
 (eta*nu)/(alpha*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
 (2*power(nu,2))/(alpha*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
 (eta*power(nu,2))/(alpha*gamma*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
 (eta*nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/
  (alpha*gamma*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
 (2*nu)/(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) + 
 (eta*nu)/(alpha*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
 (2*power(nu,2))/(alpha*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
 (eta*power(nu,2))/(alpha*gamma*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
 (eta*nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/
  (alpha*gamma*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))))/
(-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
 (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma))
 
 matrix_exp[4,2] = (exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*eta*nu*
 (-gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
 (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
(2*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
 (-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
(exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*eta*nu*
 (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
 (-gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
(2*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
 (-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
((nu*(-gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
    (gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
    (gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
  (4*alpha*power(gamma,2)*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
 (nu*(gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
    (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
    (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
  (4*alpha*power(gamma,2)*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))))/
(exp(eta*t)*(-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
   (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma))) + 
(((-gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
    (-((eta*nu*(-gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
         (alpha*gamma*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))) - 
      (nu*(gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
         (gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
       (2*alpha*gamma*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))))/(2*gamma) - 
 ((-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
    (-((eta*nu*(-gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
         (alpha*gamma*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))) - 
      (nu*(gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
         (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
       (2*alpha*gamma*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))))/(2*gamma))/
(-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
 (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma))

 matrix_exp[4,3] = 1 + (sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma - 
 (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma))/
(exp(eta*t)*(-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
   (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma)))

matrix_exp[4,4] = 1

matrix_exp[4,5] = 0
# Row 5
matrix_exp[5,1] = -(nu/(alpha - nu)) - (exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
 (-2*alpha*gamma + gamma*nu - power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
(2*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) - 
(exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
 (2*alpha*gamma - gamma*nu + power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
(2*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))
matrix_exp[5,2] = -(alpha/(alpha - nu)) + (alpha*exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
 (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
(2*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) + 
(alpha*exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
 (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
(2*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))
matrix_exp[5,3] = 0
matrix_exp[5,4] = 0
matrix_exp[5,5] = 1

return(matrix_exp)
end



# current version which works well
function newnew_eirrc_closed_solution!(outs_tmp, times, change_times, grid_size, t0, alphas, init_conds, gamma, nu, eta) 
  num_alphas = length(alphas)
  stop_times = vcat(change_times, times[end])
  start_times = vcat(times[1], change_times .+ grid_size)
  outs_matrix = get_tmp(outs_tmp, [alphas[1], gamma, nu, eta])
  current_init_conds = init_conds
  current_init_time = t0
  first_column = vcat(t0, init_conds)
  init_index = 1
  # first column is time
  # next four columns are the solutions
  for i in 1:num_alphas
    alpha = alphas[i]
    current_stop = stop_times[i]
    current_start = start_times[i]
    current_times = collect(current_start:grid_size:current_stop)
    num_times = length(current_times)
    new_indices = collect(init_index:(init_index + num_times -1))
    for j in new_indices
      new_time = current_times[j - new_indices[1] + 1]
      t = new_time - current_init_time
      @inbounds begin 
      outs_matrix[1, round(Int64, j)] = new_time

      outs_matrix[2, round(Int64, j)] = current_init_conds[1] * (-0.25*(exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
      (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
      (-2*alpha*gamma + gamma*nu - power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) - 
    (exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
     (-gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (2*alpha*gamma - gamma*nu + power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (4*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
    current_init_conds[2] * ((alpha*exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
     (-gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (4*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) + 
    (alpha*exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
     (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (4*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))

        outs_matrix[3, round(Int64, j)] = current_init_conds[1] * (-0.25*(exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
      (-2*alpha*gamma + gamma*nu - power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
      (2*alpha*gamma - gamma*nu + power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (alpha*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) - 
    (exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
     (2*alpha*gamma - gamma*nu + power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (2*alpha*gamma - gamma*nu + power(nu,2) + nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (4*alpha*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
    current_init_conds[2] * ((exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
     (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (2*alpha*gamma - gamma*nu + power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (4*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) + 
    (exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
     (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (2*alpha*gamma - gamma*nu + power(nu,2) + nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (4*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))
    

        outs_matrix[4, round(Int64, j)] = current_init_conds[1] * (-0.25*(exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*nu*
      (gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
      (gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
      (-2*alpha*gamma + gamma*nu - power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (alpha*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
      (-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
    (exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*nu*
     (gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (2*alpha*gamma - gamma*nu + power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (4*alpha*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
     (-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
    ((2*nu)/(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) - 
     (2*power(nu,2))/(alpha*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
     (2*nu)/(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) + 
     (2*power(nu,2))/(alpha*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))))/
    (exp(eta*t)*(-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
       (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma)))) + 
    current_init_conds[2] * ((exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*nu*
     (gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (4*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
     (-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
    (exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*nu*
     (gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (4*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
     (-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
    ((nu*(-gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
        (gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
        (gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
      (4*alpha*power(gamma,2)*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
     (nu*(gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
        (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
        (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
      (4*alpha*power(gamma,2)*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))))/
    (exp(eta*t)*(-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
       (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma)))) +  
    current_init_conds[3] * (-((sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma - 
      (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma))/
    (exp(eta*t)*(-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
        (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma)))))
    

      outs_matrix[5, round(Int64, j)] = current_init_conds[1] * (-0.5*(exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*eta*nu*
      (-gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
      (-2*alpha*gamma + gamma*nu - power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (alpha*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
      (-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
    (exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*eta*nu*
     (-gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (2*alpha*gamma - gamma*nu + power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (2*alpha*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
     (-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
    ((2*nu)/(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) - 
     (2*power(nu,2))/(alpha*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
     (2*nu)/(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) + 
     (2*power(nu,2))/(alpha*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))))/
    (exp(eta*t)*(-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
       (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma))) + 
    ((-2*nu)/(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) - 
     (eta*nu)/(alpha*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
     (2*power(nu,2))/(alpha*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
     (eta*power(nu,2))/(alpha*gamma*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
     (eta*nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/
      (alpha*gamma*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
     (2*nu)/(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) + 
     (eta*nu)/(alpha*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
     (2*power(nu,2))/(alpha*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
     (eta*power(nu,2))/(alpha*gamma*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
     (eta*nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/
      (alpha*gamma*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))))/
    (-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
     (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma))) + 
     current_init_conds[2] * ((exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*eta*nu*
     (-gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (2*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
     (-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
    (exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*eta*nu*
     (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (-gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (2*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
     (-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
    ((nu*(-gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
        (gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
        (gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
      (4*alpha*power(gamma,2)*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
     (nu*(gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
        (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
        (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
      (4*alpha*power(gamma,2)*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))))/
    (exp(eta*t)*(-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
       (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma))) + 
    (((-gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
        (-((eta*nu*(-gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
             (alpha*gamma*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))) - 
          (nu*(gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
             (gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
           (2*alpha*gamma*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))))/(2*gamma) - 
     ((-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
        (-((eta*nu*(-gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
             (alpha*gamma*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))) - 
          (nu*(gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
             (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
           (2*alpha*gamma*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))))/(2*gamma))/
    (-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
     (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma))) + 
    current_init_conds[3] * (1 + (sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma - 
     (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma))/
    (exp(eta*t)*(-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
       (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma)))) +
    current_init_conds[4]    

          outs_matrix[6, round(Int64, j)] = current_init_conds[1] * (-(nu/(alpha - nu)) - (exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
          (-2*alpha*gamma + gamma*nu - power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
         (2*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) - 
         (exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
          (2*alpha*gamma - gamma*nu + power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
         (2*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
         current_init_conds[2] * (-(alpha/(alpha - nu)) + (alpha*exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
          (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
         (2*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) + 
         (alpha*exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
          (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
         (2*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
         current_init_conds[5]         

      end 
    end 
    init_index = maximum(new_indices) + 1
    current_init_conds = outs_matrix[2:6,round(Int64, maximum(new_indices))]
    current_init_time = outs_matrix[1, round(Int64, maximum(new_indices))]
  end

  outs_matrix = hcat(first_column, outs_matrix)
  return(outs_matrix)
end 

function fixed_alpha_closed_soln(outs_tmp, times, grid_size, t0, alpha, init_conds, gamma, nu, eta) 
  outs_matrix = get_tmp(outs_tmp, [alphas[1], gamma, nu, eta])
  current_init_conds = init_conds
  first_column = vcat(t0, init_conds)
  # first column is time
    for j in times
      t = j - t0
      # print(t)
      @inbounds begin 
      outs_matrix[1, round(Int64, j)] = j

      outs_matrix[2, round(Int64, j)] = current_init_conds[1] * (-0.25*(exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
      (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
      (-2*alpha*gamma + gamma*nu - power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) - 
    (exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
     (-gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (2*alpha*gamma - gamma*nu + power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (4*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
    current_init_conds[2] * ((alpha*exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
     (-gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (4*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) + 
    (alpha*exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
     (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (4*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))

        outs_matrix[3, round(Int64, j)] = current_init_conds[1] * (-0.25*(exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
      (-2*alpha*gamma + gamma*nu - power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
      (2*alpha*gamma - gamma*nu + power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (alpha*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) - 
    (exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
     (2*alpha*gamma - gamma*nu + power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (2*alpha*gamma - gamma*nu + power(nu,2) + nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (4*alpha*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
    current_init_conds[2] * ((exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
     (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (2*alpha*gamma - gamma*nu + power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (4*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) + 
    (exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
     (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (2*alpha*gamma - gamma*nu + power(nu,2) + nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (4*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))
    

        outs_matrix[4, round(Int64, j)] = current_init_conds[1] * (-0.25*(exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*nu*
      (gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
      (gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
      (-2*alpha*gamma + gamma*nu - power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (alpha*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
      (-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
    (exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*nu*
     (gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (2*alpha*gamma - gamma*nu + power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (4*alpha*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
     (-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
    ((2*nu)/(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) - 
     (2*power(nu,2))/(alpha*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
     (2*nu)/(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) + 
     (2*power(nu,2))/(alpha*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))))/
    (exp(eta*t)*(-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
       (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma)))) + 
    current_init_conds[2] * ((exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*nu*
     (gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (4*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
     (-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
    (exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*nu*
     (gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (4*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
     (-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
    ((nu*(-gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
        (gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
        (gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
      (4*alpha*power(gamma,2)*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
     (nu*(gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
        (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
        (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
      (4*alpha*power(gamma,2)*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))))/
    (exp(eta*t)*(-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
       (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma)))) +  
    current_init_conds[3] * (-((sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma - 
      (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma))/
    (exp(eta*t)*(-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
        (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma)))))
    

      outs_matrix[5, round(Int64, j)] = current_init_conds[1] * (-0.5*(exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*eta*nu*
      (-gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
      (-2*alpha*gamma + gamma*nu - power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (alpha*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
      (-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
    (exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*eta*nu*
     (-gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (2*alpha*gamma - gamma*nu + power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (2*alpha*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
     (-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
    ((2*nu)/(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) - 
     (2*power(nu,2))/(alpha*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
     (2*nu)/(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) + 
     (2*power(nu,2))/(alpha*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))))/
    (exp(eta*t)*(-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
       (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma))) + 
    ((-2*nu)/(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) - 
     (eta*nu)/(alpha*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
     (2*power(nu,2))/(alpha*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
     (eta*power(nu,2))/(alpha*gamma*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
     (eta*nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/
      (alpha*gamma*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
     (2*nu)/(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) + 
     (eta*nu)/(alpha*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
     (2*power(nu,2))/(alpha*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
     (eta*power(nu,2))/(alpha*gamma*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
     (eta*nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/
      (alpha*gamma*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))))/
    (-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
     (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma))) + 
     current_init_conds[2] * ((exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*eta*nu*
     (-gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (2*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
     (-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
    (exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*eta*nu*
     (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
     (-gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
    (2*gamma*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))*
     (-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
    ((nu*(-gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
        (gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
        (gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
      (4*alpha*power(gamma,2)*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) - 
     (nu*(gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
        (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
        (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
      (4*alpha*power(gamma,2)*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))))/
    (exp(eta*t)*(-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
       (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma))) + 
    (((-gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
        (-((eta*nu*(-gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
             (alpha*gamma*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))) - 
          (nu*(gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
             (gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
           (2*alpha*gamma*(-2*eta + gamma + nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))))/(2*gamma) - 
     ((-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
        (-((eta*nu*(-gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
             (alpha*gamma*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))) - 
          (nu*(gamma - nu - sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))*
             (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
           (2*alpha*gamma*(-2*eta + gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))))/(2*gamma))/
    (-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
     (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma))) + 
    current_init_conds[3] * (1 + (sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma - 
     (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma))/
    (exp(eta*t)*(-(sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))/gamma) + 
       (nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))/(alpha*gamma)))) +
    current_init_conds[4]    

          outs_matrix[6, round(Int64, j)] = current_init_conds[1] * (-(nu/(alpha - nu)) - (exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
          (-2*alpha*gamma + gamma*nu - power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
         (2*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) - 
         (exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
          (2*alpha*gamma - gamma*nu + power(nu,2) - nu*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
         (2*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
         current_init_conds[2] * (-(alpha/(alpha - nu)) + (alpha*exp(((-gamma - sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
          (-gamma - nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
         (2*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))) + 
         (alpha*exp(((-gamma + sqrt(4*alpha*gamma + power(gamma - nu,2)) - nu)*t)/2)*
          (gamma + nu + sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2))))/
         (2*(alpha - nu)*sqrt(4*alpha*gamma + power(gamma,2) - 2*gamma*nu + power(nu,2)))) + 
         current_init_conds[5]         

      end 
    end 

  outs_matrix = hcat(first_column, outs_matrix)
  return(outs_matrix)
end 