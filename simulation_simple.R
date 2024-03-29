library(Rsolnp)
library(data.table)
library(grid)

source("hmm_functions.R")
source("hmm_parameters.R")

set.seed(321321)

params0 <- get_params0()
params1 <- get_params1(params0)
params2 <- get_params2(params0)

## This the size of the panel dataset, i.e. the number of observations per time period
n_panel_elements <- 10000

## Simulate a panel dataset
panel <- replicate(n_panel_elements, simulate_hmm(params0), simplify=FALSE)

## Each element of the panel contains an observation sequence y and hidden state sequence x
## The length of y is the number of years (or time periods) for which we observe each spatial point
panel[[1]]
length(panel[[1]]$y)

## Use the simulated panel data to estimate the model's parameters using params0, params1, and params2 as initial values

## Initialize EM at true parameter values (easy case)
params0_hat <- get_expectation_maximization_estimates(panel, params0, max_iter=20, epsilon=0.001)

stopifnot(all(diff(params0_hat$loglik) > 0))  # Loglik should be increasing
max(abs(c(params0_hat$P_list, recursive=TRUE) - c(params0$P_list, recursive=TRUE)))  # Largest error in time-varying transition probabilities
max(abs(params0_hat$pr_y - params0$pr_y))  # Largest error in observation probabilities (aka misclassification probabilities)

## Initialize EM at incorrect parameter values (more difficult)
params1_hat <- get_expectation_maximization_estimates(panel, params1, max_iter=30, epsilon=0.001)

plot(params1_hat$loglik)
max(abs(c(params1_hat$P_list, recursive=TRUE) - c(params0$P_list, recursive=TRUE)))  # Largest error in time-varying transition probabilities
max(abs(params1_hat$pr_y - params0$pr_y))  # Largest error in observation probabilities (aka misclassification probabilities)

## Initialize EM at "very" incorrect parameter values (even more difficult)
params2_hat <- get_expectation_maximization_estimates(panel, params2, max_iter=30, epsilon=0.001)

plot(params2_hat$loglik)
max(abs(c(params2_hat$P_list, recursive=TRUE) - c(params0$P_list, recursive=TRUE)))  # Largest error in time-varying transition probabilities
max(abs(params2_hat$pr_y - params0$pr_y))  # Largest error in observation probabilities (aka misclassification probabilities)

## Which initial parameter values ended up with highest likelihood? Careful, returns index in {1, 2, 3}, not {0, 1, 2}
which.max(c(max(params0_hat$loglik),
            max(params1_hat$loglik),
            max(params2_hat$loglik)))

## Convert panel from list to data.table to make it easier to compute transition matrices
for(idx in seq_along(panel)) {
    panel[[idx]]$point_id <- idx
    panel[[idx]]$time <- seq_along(panel[[idx]]$y)
}
dtable <- rbindlist(Map(data.frame, panel))
setkey(dtable, point_id)
stopifnot(all(c("point_id", "time", "x", "y") %in% names(dtable)))
table(dtable$time)  # Periods 1 through length(params0$P_list)+1

dtable[, y_one_period_ahead := c(tail(y, .N-1), NA), by="point_id"]
dtable[, y_two_periods_ahead := c(tail(y, .N-2), NA, NA), by="point_id"]
head(dtable[, c("point_id", "time", "y", "y_one_period_ahead", "y_two_periods_ahead"), with=FALSE], 25)  # Sanity check

## Compute the marginal distribution of the hidden state over time
## Compare it to the marginal distribution of observations
## "Sample" refers to the simulated panel dataset, as opposed to the population's data generating process
marginals_sample <- dtable[, list(pr_s_1=mean(x==1), pr_y_1=mean(y==1)), by="time"]
marginal_s_population <- matrix(NA, nrow=length(unique(dtable$time)), ncol=params0$n_components)
marginal_s_population[1, ] <- params0$mu
for(time_index in seq(2, nrow(marginal_s_population))) {
    marginal_s_population[time_index, ] <- marginal_s_population[time_index-1, ] %*% params0$P_list[[time_index-1]]
}
marginal_y_population <- matrix(NA, nrow=nrow(marginal_s_population), ncol=ncol(marginal_s_population))
for(time_index in seq(1, nrow(marginal_s_population))) {
    marginal_y_population[time_index, ] <- marginal_s_population[time_index, ] %*% params0$pr_y
}

## Compute the joint distribution of (Y_{t+1}, Y_{t}) for each time period t
M_Y_joint_hat_list <- lapply(seq_len(max(dtable$time) - 1), function(fixed_t) {
    with(subset(dtable, time == fixed_t), prop.table(table(y_one_period_ahead, y)))
})

## Compute inverses once and cache the results, to be re-used by solnp when estimating parameters
M_Y_joint_hat_inverse_list <- lapply(M_Y_joint_hat_list, solve)

## Compute the joint distribution of (Y_{t+2}, Y_{t+1}, Y_{t})
M_fixed_y_Y_joint_hat_list <- lapply(seq_len(params0$n_components), function(fixed_y) {
    lapply(seq_len(max(dtable$time) - 2), function(fixed_t) {
        return(with(subset(dtable, time == fixed_t & y_two_periods_ahead == fixed_y),
                    table(y_one_period_ahead, y)) / sum(dtable$time == fixed_t))
    })
})

## Naive estimate of transitions using observed Y
P1_hat_naive <- get_transition_probs_from_M_S_joint(M_Y_joint_hat_list[[1]])
P2_hat_naive <- get_transition_probs_from_M_S_joint(M_Y_joint_hat_list[[2]])
P3_hat_naive <- get_transition_probs_from_M_S_joint(M_Y_joint_hat_list[[3]])

## Compute the joint distribution of hidden states (S_t, S_{t+1}) implied by params0
M_S_joint_population <- lapply(seq_along(params0$P_list), function(time_index) {
    if(time_index == 1) {
        mu_t <- params0$mu  # Equals initial distribution when t=1
    } else {
        mu_t <- params0$mu %*% Reduce("%*%", params0$P_list[seq_len(time_index- 1)])
    }
    stopifnot(isTRUE(all.equal(sum(mu_t), 1)))  # Valid probability distribution, careful comparing floats
    return(t(params0$P_list[[time_index]] * matrix(mu_t, length(mu_t), length(mu_t))))
})
M_Y_joint_population <- lapply(M_S_joint_population, function(M) {
    return(t(params0$pr_y) %*% M %*% params0$pr_y)
})

## Compare to P1_hat_naive, P2_hat_naive, ...
P_naive_population <- lapply(M_Y_joint_population, get_transition_probs_from_M_S_joint)

M_S_joint_incorrect <- lapply(seq_along(params1$P_list), function(time_index) {
    ## Under incorrect params1 as opposed to params0
    if(time_index == 1) {
        mu_t <- params1$mu  # Equals initial distribution when t=1
    } else {
        mu_t <- params1$mu %*% Reduce("%*%", params1$P_list[seq_len(time_index- 1)])
    }
    stopifnot(isTRUE(all.equal(sum(mu_t), 1)))  # Valid probability distribution, careful comparing floats
    return(t(params1$P_list[[time_index]] * matrix(mu_t, length(mu_t), length(mu_t))))
})

## Minimum distance estimation starting from incorrect parameters
min_dist_params2_hat <- get_min_distance_estimates(params2, M_Y_joint_hat_list, M_Y_joint_hat_inverse_list, M_fixed_y_Y_joint_hat_list, dtable)
min_dist_params1_hat <- get_min_distance_estimates(params1, M_Y_joint_hat_list, M_Y_joint_hat_inverse_list, M_fixed_y_Y_joint_hat_list, dtable)

## Minimum distance estimation starting from correct parameters
min_dist_params0_hat <- get_min_distance_estimates(params0, M_Y_joint_hat_list, M_Y_joint_hat_inverse_list, M_fixed_y_Y_joint_hat_list, dtable)

## Essentially zero: we get the same min dist results starting from either params0 or params1
max(abs(min_dist_params0_hat$pr_y - min_dist_params1_hat$pr_y))
max(abs(c(min_dist_params0_hat$P_list, recursive=TRUE) - c(min_dist_params1_hat$P_list, recursive=TRUE)))

## Compare minimum distance estimates to true parameters (distances should be small)
max(abs(params0$pr_y - min_dist_params1_hat$pr_y))
max(abs(c(params0$P_list, recursive=TRUE) - c(min_dist_params1_hat$P_list, recursive=TRUE)))
max(abs(params0$mu - min_dist_params1_hat$mu))

max(abs(params0$pr_y - min_dist_params2_hat$pr_y))
max(abs(c(params0$P_list, recursive=TRUE) - c(min_dist_params2_hat$P_list, recursive=TRUE)))
max(abs(params0$mu - min_dist_params2_hat$mu))


## Check that minimum distance estimation returns correct parameter values when using population values for the distribution of Y_t
## Compare to sample analogue M_Y_joint_hat_list
## Careful, my pr_y has hidden states along the rows, transpose of misclassification matrix in paper
M_Y_joint_hat_population <- lapply(M_S_joint_population, function(M) {
    return(t(params0$pr_y) %*% M %*% params0$pr_y)
})
M_Y_joint_hat_inverse_population <- lapply(M_Y_joint_hat_population, solve)

## Compare to sample analogue M_fixed_y_Y_joint_hat_list
M_fixed_y_Y_joint_hat_population <- lapply(seq_len(params0$n_components), function(fixed_y) {
    lapply(seq_len(length(params0$P_list) - 1), function(fixed_t) {
        D <- matrix(0, params0$n_components, params0$n_components)
        diag(D) <- params0$P_list[[fixed_t + 1]] %*% t(params0$pr_y)[fixed_y, ]  # Symmetric
        return(t(params0$pr_y) %*% D %*% solve(t(params0$pr_y)) %*% M_Y_joint_hat_population[[fixed_t]])
    })
})

## Use population values for the objective function, with incorrect starting values for optimization
## Despite the incorrect starting values, the estimation should recover the true parameters without error
## The optimizer should be able to get the objective function down to zero

min_dist_params1_hat_population <- get_min_distance_estimates(params1, M_Y_joint_hat_population, M_Y_joint_hat_inverse_population, M_fixed_y_Y_joint_hat_population, dtable)

## This distance should be essentially zero: the estimates ("hats") recover the true parameters when given population (not sample) observation probabilities
max(abs(c(min_dist_params1_hat_population$P_list, recursive=TRUE) - c(params0$P_list, recursive=TRUE)))

## Optimizer gets the distance (objective function) down to zero
min_dist_params1_hat_population$objfn_values

## Try again with population, now starting the optimization from params2 ("more incorrect" than params1)
## Note convergence code of 2 and "Solution not reliable....Problem Inverting Hessian" warning
min_dist_params2_hat_population <- get_min_distance_estimates(params2, M_Y_joint_hat_population, M_Y_joint_hat_inverse_population, M_fixed_y_Y_joint_hat_population, dtable)

## Optimizer should get the distance (objective function) down to zero also when starting from params2
min_dist_params2_hat_population$objfn_values

## This distance should be zero also when starting optimization from params2
max(abs(c(min_dist_params2_hat_population$P_list, recursive=TRUE) - c(params0$P_list, recursive=TRUE)))

## Recover the initial distribution
min_dist_params2_hat_population$mu - params0$mu

## Now introduce missing observations (Y missing completely at random)
pr_missing_data <- 0.1
panel_mcar <- lapply(panel, function(panel_element) {
    panel_element_mcar <- list(time=panel_element$time,
                               point_id=panel_element$point_id,
                               x=panel_element$x)
    mask <- runif(length(panel_element$y))
    panel_element_mcar$y <- ifelse(mask < pr_missing_data, NA, panel_element$y)
    return(panel_element_mcar)
})

## Make sure EM/ML still works well with missing Ys
params0_hat_mcar <- get_expectation_maximization_estimates(panel_mcar, params0, max_iter=20, epsilon=0.001)

max(abs(c(params0_hat_mcar$P_list, recursive=TRUE) - c(params0$P_list, recursive=TRUE)))  # Largest error in time-varying transition probabilities
max(abs(params0_hat_mcar$pr_y - params0$pr_y))  # Largest error in observation probabilities (aka misclassification probabilities)

params2_hat_mcar <- get_expectation_maximization_estimates(panel_mcar, params2, max_iter=30, epsilon=0.001)

max(abs(c(params2_hat_mcar$P_list, recursive=TRUE) - c(params0$P_list, recursive=TRUE)))  # Largest error in time-varying transition probabilities
max(abs(params2_hat_mcar$pr_y - params0$pr_y))  # Largest error in observation probabilities (aka misclassification probabilities)

## Now run min dist with missing Ys
dtable_mcar <- rbindlist(Map(data.frame, panel_mcar))
setkey(dtable_mcar, point_id)
stopifnot(all(c("point_id", "time", "x", "y") %in% names(dtable_mcar)))
table(dtable_mcar$time)  # Periods 1 through length(params0$P_list)+1

dtable_mcar[, y_one_period_ahead := c(tail(y, .N-1), NA), by="point_id"]
dtable_mcar[, y_two_periods_ahead := c(tail(y, .N-2), NA, NA), by="point_id"]
head(dtable_mcar[, c("point_id", "time", "y", "y_one_period_ahead", "y_two_periods_ahead"), with=FALSE], 25)  # Sanity check

mean(is.na(dtable_mcar$y))  # This should be close to pr_missing_data

## Compute the joint distribution of (Y_{t+1}, Y_{t}) for each time period t
## Note that sum(M_Y_joint_hat_list_mcar[[1]]) and sum(M_Y_joint_hat_list[[1]]) are both 1
## prop.table(table(...)) appears to handle missing values for us (TODO Confirm!)
M_Y_joint_hat_list_mcar <- lapply(seq_len(max(dtable_mcar$time) - 1), function(fixed_t) {
    with(subset(dtable_mcar, time == fixed_t), prop.table(table(y_one_period_ahead, y)))
})

## Compute inverses once and cache the results, to be re-used by solnp when estimating parameters
M_Y_joint_hat_inverse_list_mcar <- lapply(M_Y_joint_hat_list_mcar, solve)

sum(c(M_fixed_y_Y_joint_hat_list, recursive=TRUE))  # Sum is 2.0 (because we are summing over 2 values of fixed_t)
sum(c(M_fixed_y_Y_joint_hat_list[[1]][[1]], recursive=TRUE)) + sum(c(M_fixed_y_Y_joint_hat_list[[2]][[1]], recursive=TRUE))  # Sum is 1.0

## Compute the joint distribution of (Y_{t+2}, Y_{t+1}, Y_{t})
## Note that this needs to be adjusted for MCAR relative to the simpler case when Y is never missing
M_fixed_y_Y_joint_hat_list_mcar <- lapply(seq_len(params0$n_components), function(fixed_y) {
    lapply(seq_len(max(dtable_mcar$time) - 2), function(fixed_t) {
        return(with(subset(dtable_mcar, time == fixed_t & y_two_periods_ahead == fixed_y),
                    table(y_one_period_ahead, y)) / sum(dtable_mcar$time == fixed_t &
                                                        !is.na(dtable_mcar$y_two_periods_ahead) &
                                                        !is.na(dtable_mcar$y_one_period_ahead) &
                                                        !is.na(dtable_mcar$y)))
    })
})

sum(c(M_fixed_y_Y_joint_hat_list_mcar, recursive=TRUE))  # Sum is 2.0 (because we are summing over 2 values of fixed_t)

min_dist_params2_hat_mcar <- get_min_distance_estimates(params2,
                                                        M_Y_joint_hat_list_mcar,
                                                        M_Y_joint_hat_inverse_list_mcar,
                                                        M_fixed_y_Y_joint_hat_list_mcar,
                                                        dtable_mcar)

max(abs(params0$pr_y - min_dist_params2_hat_mcar$pr_y))
max(abs(c(params0$P_list, recursive=TRUE) - c(min_dist_params2_hat_mcar$P_list, recursive=TRUE)))
max(abs(params0$mu - min_dist_params2_hat_mcar$mu))
