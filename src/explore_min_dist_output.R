max_abs <- function(x) {
    return(max(abs(x)))
}

## Note: this file was generated _before_ adding diag dominant pr_y logic to hmm_functions.R
infile <- "county_simulation_2019-10-15 20:10:08_iter_3.rds"

## Note: this file was generated _after_ adding diag dominant pr_y logic to hmm_functions.R
infile <- "county_simulation_2019-10-16 18:51:25_iter_2.rds"

simulations <- readRDS(infile)

## Difference between min dist estimates of P and true P
P_errors <- lapply(simulations, function(x) {
    c(x$estimates$min_dist_params_hat_best_objfn$P_list, recursive=TRUE) - c(x$params$P_list, recursive=TRUE)
})

hist(sapply(P_errors, mean))  # Looks good, means are close to zero
hist(sapply(P_errors, sd))  # A few std devs are worryingly large
hist(sapply(P_errors, max_abs))  # Some large errors in min dist estimates of P

max_abs_errors <- sapply(P_errors, max_abs)

## Note: this fraction is much larger with n_points_per_county=100 than with n_points_per_county=500
mean(max_abs_errors > 0.5)
indices_large_errors <- which(max_abs_errors > 0.5)

simulations_large_errors <- simulations[indices_large_errors]

simulations_large_errors[[1]]$params

simulations_large_errors[[1]]$estimates$min_dist_params_hat_best_objfn

## Something is wrong: this isn't a valid pr_y matrix (it's not diagonally dominant)
## Compare to simulations_large_errors[[1]]$params$pr_y
simulations_large_errors[[1]]$estimates$min_dist_params_hat_best_objfn$pr_y

## Compare to simulations_large_errors[[1]]$params$mu
simulations_large_errors[[1]]$estimates$min_dist_params_hat_best_objfn$mu

## EM did fine in this example
simulations_large_errors[[1]]$estimates$em_params_hat_best_likelihood$pr_y
simulations_large_errors[[1]]$estimates$em_params_hat_best_likelihood$P_list

## Note: this is NULL when all min dist estimates had pr_y that are not diagonally dominant
simulations_large_errors[[1]]$estimates$min_dist_params_hat_diag_dominant_best_objfn
