max_abs <- function(x) {
    return(max(abs(x)))
}

infile <- "county_simulation_2019-10-15 20:10:08_iter_3.rds"

simulations <- readRDS(infile)

## Difference between min dist estimates of P and true P
P_errors <- lapply(simulations, function(x) {
    c(x$estimates$min_dist_params_hat_best_objfn$P_list, recursive=TRUE) - c(x$params$P_list, recursive=TRUE)
})

hist(sapply(P_errors, mean))  # Looks good, means are close to zero
hist(sapply(P_errors, sd))  # A few std devs are worryingly large
hist(sapply(P_errors, max_abs))  # Some large errors in min dist estimates of P

max_abs_errors <- sapply(P_errors, max_abs)
indices_large_errors <- which(max_abs_errors > 0.5)

simulations_large_errors <- simulations[indices_large_errors]

simulations_large_errors[[1]]$params

simulations_large_errors[[1]]$estimates$min_dist_params_hat_best_objfn

## Something is wrong: this isn't a valid pr_y matrix (it's not diagonally dominant)
simulations_large_errors[[1]]$estimates$min_dist_params_hat_best_objfn$pr_y

## This also looks wrong: mu is a scalar but should be a vector
## Compare to simulations_large_errors[[1]]$params$mu
simulations_large_errors[[1]]$estimates$min_dist_params_hat_best_objfn$mu

## EM did fine in this example
simulations_large_errors[[1]]$estimates$em_params_hat_best_likelihood$pr_y
