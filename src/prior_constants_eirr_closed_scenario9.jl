# EIRR-ww prior parameters
# scenario 9, misspecify E and I, and R1 to be centered at 0.75% of true values
const gamma_sd = 0.2
const gamma_mean =log(1/4)
const nu_sd = 0.2
const nu_mean = log(1/7)
const eta_sd = 0.2
const eta_mean = log(1/18)
const rho_gene_sd = 1
const rho_gene_mean = 0
const tau_sd = 1
const tau_mean = 0
const I_init_sd = 0.05
const I_init_mean = 489 * 0.75
const R1_init_sd = 0.05
const R1_init_mean = 2075 * 0.75
const E_init_sd = 0.05
const E_init_mean = 225 * 0.75
const lambda_mean = 5.685528
const lambda_sd = 2.178852
const df_shape = 2
const df_scale = 10
const sigma_rt_sd = 0.2
const sigma_rt_mean = log(0.1)
const rt_init_sd = 0.1
const rt_init_mean = log(0.88)
