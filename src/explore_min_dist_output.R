source("hmm_functions.R")

max_abs <- function(x) {
    return(max(abs(x)))
}

## Note: this file was generated _before_ adding diag dominant pr_y constraint to hmm_functions.R
infile <- "county_simulation_2019-10-15 20:10:08_iter_3.rds"

## Generated _after_ adding diag dominant pr_y constraint to hmm_functions.R
infile <- "county_simulation_2019-10-17 18:47:20_iter_1.rds"  # Only 100 points per county (tiny sample size!)
infile <- "county_simulation_2019-10-17 18:47:20_iter_2.rds"  # 500 points per county

simulations <- readRDS(infile)

## Difference between min dist estimates of P and true P
P_errors_min_dist <- lapply(simulations, function(x) {
    sapply(x$estimates$min_dist_params_hat_best_objfn$P_list, get_deforestation_prob_from_P) - sapply(x$params$P_list, get_deforestation_prob_from_P)
})

hist(sapply(P_errors_min_dist, mean))  # Biased upwards -- makes sense if estimates hit edge of parameter space with some non-negligible probability
mean(c(P_errors_min_dist, recursive=TRUE))  # 0.036
sd(c(P_errors_min_dist, recursive=TRUE))  # 0.11
hist(sapply(P_errors_min_dist, max_abs))  # Some large errors in min dist estimates of P

P_errors_em <- lapply(simulations, function(x) {
    sapply(x$estimates$em_params_hat_best_likelihood$P_list, get_deforestation_prob_from_P) - sapply(x$params$P_list, get_deforestation_prob_from_P)
})

hist(sapply(P_errors_em, mean))  # Also biased upwards
mean(c(P_errors_em, recursive=TRUE))  # 0.02, smaller than min dist
sd(c(P_errors_em, recursive=TRUE))  # 0.07
hist(sapply(P_errors_em, max_abs))  # Some large errors in em estimates of P, but better than min dist

## For each simulation, what was the largest error in min dist estimates of the deforestation rate?
max_abs_errors <- sapply(P_errors_min_dist, max_abs)

## Note: this fraction is much larger with n_points_per_county=100 than with n_points_per_county=500
mean(max_abs_errors > 0.2)
indices_large_errors <- which(max_abs_errors > 0.2)

simulations_large_errors <- simulations[indices_large_errors]

simulations_large_errors[[1]]$params

simulations_large_errors[[1]]$estimates$min_dist_params_hat_best_objfn

## Compare to simulations_large_errors[[1]]$params$pr_y
## Note: the min dist estimate of pr_y can hit the diag dominant constraint (more likely at small sample sizes?)
simulations_large_errors[[1]]$estimates$min_dist_params_hat_best_objfn$pr_y

## Compare to simulations_large_errors[[1]]$params$mu
simulations_large_errors[[1]]$estimates$min_dist_params_hat_best_objfn$mu

## EM did fine (or at least better than min dist) in this example
simulations_large_errors[[1]]$estimates$em_params_hat_best_likelihood$pr_y
simulations_large_errors[[1]]$estimates$em_params_hat_best_likelihood$P_list
