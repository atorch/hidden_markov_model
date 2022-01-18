## TODO Run with more simulations (50, 100?)
## Then check RMSEs

## TODO Make sure I'm saving MD runtimes, need to account for those (in MD+ML)

## Value used in simulation_spatial_corr.R for P_list[[1]][1, 2]
true_deforestation_rate <- 0.04

## EM using MD as initial value
sims = readRDS("spatial_corr_n_fields_10000_pr_missing_data_0_include_z_FALSE_ising_beta_0_25_simulations_use_md_as_initial_values_for_em_TRUE.rds")
em_runtimes_using_md_as_initial_values <- sapply(sims, function(sim) {
    em_runtime <- (as.POSIXct(sim$estimates$em_params_hat_list[[1]]$time_finished_em) -
                   as.POSIXct(sim$estimates$em_params_hat_list[[1]]$time_started_em))
    em_runtime <- as.numeric(em_runtime, units="secs")
    return(em_runtime)
})
em_iterations_using_md_as_initial_values <- sapply(sims, function(sim) {
    return(sim$estimates$em_params_hat_list[[1]]$n_em_iterations)
})
em_squared_error_using_md_as_initial_values <- sapply(sims, function(sim) {
    estimated_deforestation_rate <- sim$estimates$em_params_hat_list[[1]]$P_list[[1]][1, 2]
    return((estimated_deforestation_rate - true_deforestation_rate)^2)
})
md_runtimes <- sapply(sims, function(sim) {
    md_runtime <- (as.POSIXct(sim$estimates$min_dist_params_hat[[1]]$time_finished_md) -
                   as.POSIXct(sim$estimates$min_dist_params_hat[[1]]$time_started_md))
    md_runtime <- as.numeric(md_runtime, units="secs")
    return(md_runtime)
})
md_squared_error <- sapply(sims, function(sim) {
    estimated_deforestation_rate <- sim$estimates$min_dist_params_hat[[1]]$P_list[[1]][1, 2]
    return((estimated_deforestation_rate - true_deforestation_rate)^2)
})

## EM using random initial values
sims = readRDS("spatial_corr_n_fields_10000_pr_missing_data_0_include_z_FALSE_ising_beta_0_25_simulations_use_md_as_initial_values_for_em_FALSE.rds")
em_runtimes_using_random_initial_values <- sapply(sims, function(sim) {
    em_runtime <- (as.POSIXct(sim$estimates$em_params_hat_list[[1]]$time_finished_em) -
                   as.POSIXct(sim$estimates$em_params_hat_list[[1]]$time_started_em))
    em_runtime <- as.numeric(em_runtime, units="secs")
    return(em_runtime)
})
em_iterations_using_random_initial_values <- sapply(sims, function(sim) {
    return(sim$estimates$em_params_hat_list[[1]]$n_em_iterations)    
})
em_squared_error_using_random_initial_values <- sapply(sims, function(sim) {
    estimated_deforestation_rate <- sim$estimates$em_params_hat_list[[1]]$P_list[[1]][1, 2]
    return((estimated_deforestation_rate - true_deforestation_rate)^2)
})

t.test(em_runtimes_using_random_initial_values, md_runtimes)
t.test(em_runtimes_using_random_initial_values, md_runtimes + em_runtimes_using_md_as_initial_values)

t.test(em_squared_error_using_random_initial_values, md_squared_error)

t.test(em_runtimes_using_md_as_initial_values, em_runtimes_using_random_initial_values)
t.test(em_iterations_using_md_as_initial_values, em_iterations_using_random_initial_values)
t.test(em_squared_error_using_md_as_initial_values, em_squared_error_using_random_initial_values)
