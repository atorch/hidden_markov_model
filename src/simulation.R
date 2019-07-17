library(data.table)
library(ggplot2)
library(grid)
library(latex2exp)  # For ggplot2 xlab
library(parallel)
library(Rsolnp)

## Variables defined in this simulation can be reference in latex documents using \Sexpr

source("functions_hmm.R")
source("ggplot_utils.R")
set_ggplot_theme()

set.seed(321321)

get_random_initial_parameters <- function(params0) {

    ## Given a true set of HMM parameters, return random incorrect parameters from which to begin parameter estimation

    initial_parameters <- list(n_components=params0$n_components)

    initial_parameters$P_list <- lapply(params0$P_list, function(correct_P) {

        ## Probabilities on diagonals of the transition probability matrices
        random_uniform <- runif(params0$n_components, min=0.60, max=0.98)

        P <- matrix((1 - random_uniform) / (params0$n_components - 1), nrow=nrow(correct_P), ncol=ncol(correct_P))
        diag(P) <- random_uniform

        return(P)
    })

    ## Probabilities on the diagonals of the observation probability matrix pr_y
    random_uniform <- runif(params0$n_components, min=0.60, max=0.98)
    initial_parameters$pr_y <- matrix((1 - random_uniform) / (params0$n_components - 1), nrow=nrow(params0$pr_y), ncol=ncol(params0$pr_y))
    diag(initial_parameters$pr_y) <- random_uniform

    ## The initial distribution over hidden states is set to its true value
    initial_parameters$mu <- params0$mu

    return(initial_parameters)
}

get_hmm_and_minimum_distance_estimates_random_initialization <- function(params0, n_panel_elements=5000, n_random_starts=5) {

    ## Params0 are true parameters, params1 are incorrect parameters from which to start estimation

    require(data.table)
    require(Rsolnp)

    panel <- replicate(n_panel_elements, simulate_hmm(params0), simplify=FALSE)
    random_initial_parameters <- replicate(n=n_random_starts, get_random_initial_parameters(params0), simplify=FALSE)

    hmm_params_hat_list <- lapply(random_initial_parameters, function(initial_params) {
        return(em_parameter_estimates(panel, initial_params, max_iter=30, epsilon=0.001))
    })
    likelihoods <- sapply(hmm_params_hat_list, function(x) {
        return(max(x$loglik))
    })

    for(idx in seq_along(panel)) {
        panel[[idx]]$point_id <- idx
        panel[[idx]]$time <- seq_along(panel[[idx]]$y)
    }
    dtable <- rbindlist(Map(data.frame, panel))
    setkey(dtable, point_id)
    stopifnot(all(c("point_id", "time", "x", "y") %in% names(dtable)))
    dtable[, y_one_period_ahead := c(tail(y, .N-1), NA), by="point_id"]
    dtable[, y_two_periods_ahead := c(tail(y, .N-2), NA, NA), by="point_id"]
    head(dtable[, c("point_id", "time", "y", "y_one_period_ahead", "y_two_periods_ahead"), with=FALSE], 25)  # Sanity check
    M_Y_joint_hat_list <- lapply(seq_len(max(dtable$time) - 1), function(fixed_t) {
        with(subset(dtable, time == fixed_t), prop.table(table(y_one_period_ahead, y)))
    })  # Joint distribution of (Y_{t+1}, Y_{t}) and (Y_{t+2}, Y_{t+1})
    M_Y_joint_hat_inverse_list <- lapply(M_Y_joint_hat_list, solve)  # Compute inverses once, before running solnp
    M_fixed_y_Y_joint_hat_list <- lapply(seq_len(params0$n_components), function(fixed_y) {
        lapply(seq_len(max(dtable$time) - 2), function(fixed_t) {
            return(with(subset(dtable, time == fixed_t & y_two_periods_ahead == fixed_y),
                        table(y_one_period_ahead, y)) / sum(dtable$time == fixed_t))
        })
    })  # TODO Need to handle edge case where any of these matrices are not invertible, might happen at small sample sizes
    min_dist_params_hat_list <- lapply(random_initial_parameters, function(initial_params) {
        M_S_joint_list_initial <- lapply(seq_along(initial_params$P_list), function(time_index) {
            ## Joint distribution of S_t, S_{t+1} implied by initial params
            if(time_index == 1) {
                mu_t <- initial_params$mu  # Equals initial distribution when t=1
            } else {
                mu_t <- initial_params$mu %*% Reduce("%*%", initial_params$P_list[seq_len(time_index- 1)])
            }
            stopifnot(isTRUE(all.equal(sum(mu_t), 1)))  # Valid probability distribution, careful comparing floats
            return(t(initial_params$P_list[[time_index]] * matrix(mu_t, length(mu_t), length(mu_t))))
        })
        x_guess1 <- c(t(initial_params$pr_y), c(M_S_joint_list_initial, recursive=TRUE))
        max_time <- max(dtable$time)
        solnp_result1 <- solnp(x_guess1,
                               fun=objfn_minimum_distance, eqfun=eq_function_minimum_distance,
                               eqB=rep(1, params0$n_components + max(dtable$time) - 1),
                               LB=rep(0, length(x_guess1)),
                               UB=rep(1, length(x_guess1)),
                               M_Y_joint_hat_list=M_Y_joint_hat_list,
                               M_Y_joint_hat_inverse_list=M_Y_joint_hat_inverse_list,
                               M_fixed_y_Y_joint_hat_list=M_fixed_y_Y_joint_hat_list,
                               max_time=max_time,
                               n_components=params0$n_components,
                               control=list(delta=1e-14, tol=1e-14, trace=1))  # Careful, sensitive to control
        M_Y_given_S_hat1 <- matrix(solnp_result1$pars[seq(1, params0$n_components^2)], params0$n_components, params0$n_components)  # Transpose of params0$pr_y
        M_S_joint_list_hat1 <- lapply(seq_len(max_time - 1), function(time_index, n_components=initial_params$n_components) {
            return(matrix(solnp_result1$pars[seq((n_components^2)*time_index + 1, (n_components^2)*(1 + time_index))], n_components, n_components))

        })
        min_dist_params_hat <- list(pr_y=t(M_Y_given_S_hat1),
                                    P_list=lapply(M_S_joint_list_hat1, get_transition_probs_from_M_S_joint),
                                    convergence=solnp_result1$convergence,
                                    objfn_values=solnp_result1$values)  # Choose estimates with lowest objfn values
        return(min_dist_params_hat)
    })
    objfn_values <- sapply(min_dist_params_hat_list, function(x) {
        return(min(x$objfn_values))
    })
    return(list("panel_size"=n_panel_elements,
                "hmm_params_hat_list"=hmm_params_hat_list,
                "hmm_params_hat_loglikelihoods"=likelihoods,
                "initial_parameters_list"=random_initial_parameters,
                "hmm_params_hat_best_likelihood"=hmm_params_hat_list[[which.max(likelihoods)]],
                "min_dist_params_hat_list"=min_dist_params_hat_list,
                "min_dist_objfn_values"=objfn_values,
                "min_dist_params_hat_best_objfn"=min_dist_params_hat_list[[which.min(objfn_values)]]))
}

get_params0 <- function() {

    ## Build a list of parameters defining an HMM model, used as input in simulation and estimation functions

    ## Initial distribution over hidden states
    params0 <- list(n_components=2, mu=c(0.7, 0.3))

    ## Observation probabilities (rows are hidden states, columns are Y)
    pr_y <- rbind(c(0.90, 0.10),
                  c(0.20, 0.80))

    ## List of transition probabilities for the hidden state
    ## P_list[[t]] gives the transition probabilties from period t to t+1
    P_list <- list(rbind(c(0.96, 0.04),
                         c(0.02, 0.98)),
                   rbind(c(0.90, 0.10),
                         c(0.07, 0.93)),
                   rbind(c(0.80, 0.20),
                         c(0.30, 0.70)))

    stopifnot(rows_sum_to_one(pr_y))
    stopifnot(all(sapply(P_list, rows_sum_to_one)))

    params0$P_list <- P_list
    params0$pr_y <- pr_y

    return(params0)
}

get_params1 <- function(params0) {

    ## Define incorrect parameters, see whether HMM and min distance estimators still do well under incorrect starting values
    params1 <- params0

    ## Modify the observation probabilities
    params1$pr_y <- rbind(c(0.70, 0.30),
                          c(0.30, 0.70))

    ## Modify the transition probabilities
    params1$P_list <- replicate(n=length(params0$P_list),
                                rbind(c(0.96, 0.04),
                                      c(0.04, 0.96)), simplify=FALSE)

    return(params1)
}

get_params2 <- function(params0) {

    ## Define second set of incorrect parameters, even farther from params0 than params1
    params2 <- params0

    ## Modify the observation probabilities
    params2$pr_y <- rbind(c(0.60, 0.40),
                          c(0.40, 0.60))

    ## Modify the transition probabilities
    params2$P_list <- replicate(n=length(params0$P_list),
                                rbind(c(0.60, 0.40),
                                      c(0.40, 0.60)), simplify=FALSE)

    return(params2)
}

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

outfile <- "simulation_methodology_params0_hat.rds"
if(file.exists(outfile)) {
    params0_hat <- readRDS(outfile)
} else {
    params0_hat <- em_parameter_estimates(panel, params0, max_iter=20, epsilon=0.001)  # Starting from true parameter values
    saveRDS(params0_hat, file=outfile)
}
stopifnot(all(diff(params0_hat$loglik) > 0))  # Loglik should be increasing
max(abs(c(params0_hat$P_list, recursive=TRUE) - c(params0$P_list, recursive=TRUE)))  # Largest error in time-varying transition probabilities
max(abs(params0_hat$pr_y - params0$pr_y))  # Largest error in observation probabilities (aka misclassification probabilities)

outfile <- "simulation_methodology_params1_hat.rds"
if(file.exists(outfile)) {
    params1_hat <- readRDS(outfile)
} else {
    params1_hat <- em_parameter_estimates(panel, params1, max_iter=20, epsilon=0.001)  # Starting from incorrect parameter values
    saveRDS(params1_hat, file=outfile)
}
max(abs(c(params1_hat$P_list, recursive=TRUE) - c(params0$P_list, recursive=TRUE)))  # Largest error in time-varying transition probabilities
max(abs(params1_hat$pr_y - params0$pr_y))  # Largest error in observation probabilities (aka misclassification probabilities)

outfile <- "simulation_methodology_params2_hat.rds"
if(file.exists(outfile)) {
    params2_hat <- readRDS(outfile)
} else {
    params2_hat <- em_parameter_estimates(panel, params2, max_iter=20, epsilon=0.001)  # Starting from "very" incorrect parameter values
    saveRDS(params2_hat, file=outfile)
}
max(abs(c(params2_hat$P_list, recursive=TRUE) - c(params0$P_list, recursive=TRUE)))  # Largest error in time-varying transition probabilities
max(abs(params2_hat$pr_y - params0$pr_y))  # Largest error in observation probabilities (aka misclassification probabilities)

## Which initial parameter values ended up with highest likelihood (careful, returns index in {1, 2, 3}, not {0, 1, 2})
which.max(c(max(params0_hat$loglik),
            max(params1_hat$loglik),
            max(params2_hat$loglik)))

## Convert panel from list to data.table, to make it easier to compute transition matrices
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

## Compute the joint distribution of (Y_{t+1}, Y_{t}) conditional on Y_{t+2}
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
M_S_joint_list_population <- lapply(seq_along(params0$P_list), function(time_index) {
    if(time_index == 1) {
        mu_t <- params0$mu  # Equals initial distribution when t=1
    } else {
        mu_t <- params0$mu %*% Reduce("%*%", params0$P_list[seq_len(time_index- 1)])
    }
    stopifnot(isTRUE(all.equal(sum(mu_t), 1)))  # Valid probability distribution, careful comparing floats
    return(t(params0$P_list[[time_index]] * matrix(mu_t, length(mu_t), length(mu_t))))
})
M_Y_joint_list_population <- lapply(M_S_joint_list_population, function(M) {
    return(t(params0$pr_y) %*% M %*% params0$pr_y)
})

## Compare to P1_hat_naive, P2_hat_naive, ...
P_naive_population <- lapply(M_Y_joint_list_population, get_transition_probs_from_M_S_joint)

## Minimum distance estimation starting from incorrect parameters
M_S_joint_list_incorrect <- lapply(seq_along(params1$P_list), function(time_index) {
    ## Under incorrect params1 as opposed to params0
    if(time_index == 1) {
        mu_t <- params1$mu  # Equals initial distribution when t=1
    } else {
        mu_t <- params1$mu %*% Reduce("%*%", params1$P_list[seq_len(time_index- 1)])
    }
    stopifnot(isTRUE(all.equal(sum(mu_t), 1)))  # Valid probability distribution, careful comparing floats
    return(t(params1$P_list[[time_index]] * matrix(mu_t, length(mu_t), length(mu_t))))
})

x_guess1 <- c(t(params1$pr_y), c(M_S_joint_list_incorrect, recursive=TRUE))  # Incorrect HMM parameters as initial guess
max_time <- max(dtable$time)
solnp_result1 <- solnp(x_guess1,
                       fun=objfn_minimum_distance, eqfun=eq_function_minimum_distance,
                       eqB=rep(1, params0$n_components + max(dtable$time) - 1),
                       LB=rep(0, length(x_guess1)),
                       UB=rep(1, length(x_guess1)),
                       M_Y_joint_hat_list=M_Y_joint_hat_list,
                       M_Y_joint_hat_inverse_list=M_Y_joint_hat_inverse_list,
                       M_fixed_y_Y_joint_hat_list=M_fixed_y_Y_joint_hat_list,
                       max_time=max_time,
                       n_components=params0$n_components,
                       control=list(delta=1e-14, tol=1e-14, trace=1))  # Careful, sensitive to control
M_Y_given_S_hat1 <- matrix(solnp_result1$pars[seq(1, params0$n_components^2)], params0$n_components, params0$n_components)  # Transpose of params0$pr_y
M_S_joint_list_hat1 <- lapply(seq_len(max_time - 1), function(time_index, n_components=params1$n_components) {
    return(matrix(solnp_result1$pars[seq((n_components^2)*time_index + 1, (n_components^2)*(1 + time_index))], n_components, n_components))
})
min_dist_params1_hat <- list(pr_y=t(M_Y_given_S_hat1),
                             P_list=lapply(M_S_joint_list_hat1, get_transition_probs_from_M_S_joint),
                             convergence=solnp_result1$convergence,
                             objfn_values=solnp_result1$values)

## Minimum distance estimation starting from correct parameters
x_guess0 <- c(t(params0$pr_y), c(M_S_joint_list_population, recursive=TRUE))
solnp_result0 <- solnp(x_guess0,
                       fun=objfn_minimum_distance, eqfun=eq_function_minimum_distance,
                       eqB=rep(1, params0$n_components + max(dtable$time) - 1),
                       LB=rep(0, length(x_guess1)),
                       UB=rep(1, length(x_guess1)),
                       M_Y_joint_hat_list=M_Y_joint_hat_list,
                       M_Y_joint_hat_inverse_list=M_Y_joint_hat_inverse_list,
                       M_fixed_y_Y_joint_hat_list=M_fixed_y_Y_joint_hat_list,
                       max_time=max_time,
                       n_components=params0$n_components,
                       control=list(delta=1e-14, tol=1e-14, trace=1))  # Careful, sensitive to control
M_Y_given_S_hat0 <- matrix(solnp_result0$pars[seq(1, params0$n_components^2)], params0$n_components, params0$n_components)  # Transpose of params0$pr_y
M_S_joint_list_hat0 <- lapply(seq_len(max_time - 1), function(time_index, n_components=params0$n_components) {
    return(matrix(solnp_result0$pars[seq((n_components^2)*time_index + 1, (n_components^2)*(1 + time_index))], n_components, n_components))
})
min_dist_params0_hat <- list(pr_y=t(M_Y_given_S_hat0),
                             P_list=lapply(M_S_joint_list_hat0, get_transition_probs_from_M_S_joint),
                             convergence=solnp_result0$convergence,
                             objfn_values=solnp_result0$values)

## Essentially zero: we get the same min dist results starting from either params0 or params1
max(abs(M_Y_given_S_hat0 - M_Y_given_S_hat1))
max(abs(c(M_S_joint_list_hat0, recursive=TRUE) - c(M_S_joint_list_hat1, recursive=TRUE)))

## Check that minimum distance estimation returns correct parameter values when using population values for the distribution of Y_t
## Compare to sample analogue M_Y_joint_hat_list
M_Y_joint_hat_list_population <- lapply(M_S_joint_list_population, function(M) {
    return(t(params0$pr_y) %*% M %*% params0$pr_y)  # Careful, my pr_y has hidden states along the rows, transpose of misclassification matrix in paper
})
M_Y_joint_hat_inverse_list_population <- lapply(M_Y_joint_hat_list_population, solve)
M_fixed_y_Y_joint_hat_list_population <- lapply(seq_len(params0$n_components), function(fixed_y) {
    lapply(seq_len(length(params0$P_list) - 1), function(fixed_t) {
        D <- matrix(0, params0$n_components, params0$n_components)
        diag(D) <- params0$P_list[[fixed_t + 1]] %*% t(params0$pr_y)[fixed_y, ]  # Symmetric
        return(t(params0$pr_y) %*% D %*% solve(t(params0$pr_y)) %*% M_Y_joint_hat_list_population[[fixed_t]])  # Compare to sample analogue M_fixed_y_Y_joint_hat_list
    })
})
## Use population values for the objective function, with incorrect starting values for optimization
## Despite the incorrect starting values, the estimation should recover the true parameters without error
## The optimizer should be able to get the objective function down to zero
solnp_result_population <- solnp(x_guess1,
                                 fun=objfn_minimum_distance, eqfun=eq_function_minimum_distance,
                                 eqB=rep(1, params0$n_components + max(dtable$time) - 1),
                                 LB=rep(0, length(x_guess1)),
                                 UB=rep(1, length(x_guess1)),
                                 M_Y_joint_hat_list=M_Y_joint_hat_list_population,
                                 M_Y_joint_hat_inverse_list=M_Y_joint_hat_inverse_list_population,
                                 M_fixed_y_Y_joint_hat_list=M_fixed_y_Y_joint_hat_list_population,
                                 max_time=max_time,
                                 n_components=params0$n_components,
                                 control=list(delta=1e-14, tol=1e-14, trace=1)) 

## Compare to true parameters i.e. t(params0$pr_y)
M_Y_given_S_hat_population <- matrix(solnp_result_population$pars[seq(1, params0$n_components^2)], params0$n_components, params0$n_components)

## Compare to M_S_joint_list_population
M_S_joint_list_hat_population <- lapply(seq_len(max_time - 1), function(time_index, n_components=params1$n_components) {
    return(matrix(solnp_result_population$pars[seq((n_components^2)*time_index + 1, (n_components^2)*(1 + time_index))], n_components, n_components))
})

## This distance should be essentially zero: the estimates ("hats") recover the true parameters
max(abs(c(M_S_joint_list_hat_population, recursive=TRUE) - c(M_S_joint_list_population, recursive=TRUE)))

## Second simulation using random initial parameter values
test_random_params <- get_random_initial_parameters(params0)
## test_estimates_with_5_random_initializations <- get_hmm_and_minimum_distance_estimates_random_initialization(params0, n_random_starts=5)  # Slow!

cluster <- makeCluster(detectCores())
n_replications <- 100
n_random_starts <- 6
outfile_format <- "simulation_sampling_distribution_with_%s_random_initial_parameters_panel_size_%s_%s_replications.rds"
clusterExport(cluster, c("get_random_initial_parameters",
                         "get_hmm_and_minimum_distance_estimates_random_initialization",
                         "em_parameter_estimates",
                         "baum_welch",
                         "objfn_minimum_distance",
                         "eq_function_minimum_distance",
                         "get_transition_probs_from_M_S_joint",
                         "simulate_hmm",
                         "simulate_discrete_markov",
                         "valid_panel_element",
                         "valid_parameters",
                         "params0",
                         "params1",
                         "n_random_starts"))

panel_sizes <- c(100, 200, 500, 1000, 5000, 10000)  # Slow, panel size 5000 took 22 hours on my laptop (4 CPUs)
panel_sizes <- c(100, 200, 500)

for(panel_size in panel_sizes) {
    outfile <- sprintf(outfile_format, n_random_starts, panel_size, n_replications)
    if(file.exists(outfile)) {
        message(outfile, " already exists; won't re-run simulation")
    } else {
        message("running replications with random initialization, time is ", Sys.time())
        replications_random_initialization <- parLapply(cluster, rep(panel_size, n_replications), function(x) {
            get_hmm_and_minimum_distance_estimates_random_initialization(params0=params0, n_panel_elements=x, n_random_starts=n_random_starts)  # Slow...
        })
        message("done, time is ", Sys.time())
        message("saving ", outfile)
        saveRDS(replications_random_initialization, outfile)
    }
}
stopCluster(cluster)

infiles <- sprintf(outfile_format, n_random_starts, panel_sizes, n_replications)

## This is a list of lists: for each set of parameter values, we have multiple replications of the simulation
all_simulations <- lapply(infiles, readRDS)
stopifnot(all(sapply(all_simulations, length) == n_replications))

dataframes <- lapply(all_simulations, function(replications) {
    df <- data.frame(panel_size=sapply(replications, function(replication) {
        return(replication$hmm_params_hat_best_likelihood$panel_size)
    }),
    hmm_params_hat_n_iterations=sapply(replications, function(replication) {
        return(replication$hmm_params_hat_best_likelihood$n_em_iterations)
    }),
    hmm_params_hat_pr_y_11=sapply(replications, function(replication) {
        return(replication$hmm_params_hat_best_likelihood$pr_y[1, 1])
    }),
    hmm_params_hat_pr_y_22=sapply(replications, function(replication) {
        return(replication$hmm_params_hat_best_likelihood$pr_y[2, 2])
    }),
    min_dist_params_hat_pr_y_11=sapply(replications, function(replication) {
        return(replication$min_dist_params_hat_best_objfn$pr_y[1, 1])
    }),
    min_dist_params_hat_pr_y_22=sapply(replications, function(replication) {
        return(replication$min_dist_params_hat_best_objfn$pr_y[2, 2])
    }))
    for(time_index in seq_along(params0$P_list)) {
        varname <- sprintf("hmm_params_hat_P%s_11", time_index)
        df[, varname] <- sapply(replications, function(replication) {
            return(replication$hmm_params_hat_best_likelihood$P_list[[time_index]][1, 1])
        })
        varname <- sprintf("min_dist_params_hat_P%s_11", time_index)
        df[, varname] <- sapply(replications, function(replication) {
            return(replication$min_dist_params_hat_best_objfn$P_list[[time_index]][1, 1])
        })
    }
    return(df)
})

df <- rbindlist(dataframes)
df$panel_size_label <- sprintf("%s%s", ifelse(df$panel_size %in% c(100, 1000), "sample size = ", ""), df$panel_size)
df$panel_size_label <- factor(df$panel_size_label,
                              levels=sprintf("%s%s", ifelse(sort(unique(df$panel_size)) %in% c(100, 1000), "sample size = ", ""), sort(unique(df$panel_size))))

## Transition probs
for (time_index in seq_along(params0$P_list)) {
    for(estimator in c("hmm", "min_dist")) {  # HMM means EM
        variable <- sprintf("%s_params_hat_P%s_11", estimator, time_index)
        p <- (ggplot(df, aes_string(x=variable)) +
              geom_vline(xintercept=params0$P_list[[time_index]][1, 1], color="grey", linetype=2) +
              geom_histogram(binwidth=0.005, color="black", fill="white") +
              facet_wrap(~ panel_size_label, scale="free_x") +
              geom_rug(aes(x=true_value), data=data.frame(true_value=params0$P_list[[time_index]][1, 1])) +  # Correct parameter value
              ylab("") + xlab(""))
        outfile <- sprintf("simulation_methodology_%s_random_initialization_P_hat_%s_11_sampling_distribution.png", estimator, time_index)
        message("saving ", outfile)
        ggsave(outfile, p, width=12.5, height=8)
    }
}

## Misclassification probs
for(matrix_index in c("11", "22")) {
    for(estimator in c("hmm", "min_dist")) {  # HMM means EM
        variable <- sprintf("%s_params_hat_pr_y_%s", estimator, matrix_index)
        true_value <- ifelse(matrix_index == "11", params0$pr_y[1, 1], params0$pr_y[2, 2])
        p <- (ggplot(df, aes_string(x=variable)) +
              geom_vline(xintercept=true_value, color="grey", linetype=2) +
              geom_histogram(binwidth=0.005, color="black", fill="white") +
              facet_wrap(~ panel_size_label, scale="free_x") +
              geom_rug(aes(x=true_value), data=data.frame(true_value=true_value)) +  # Correct parameter value
              ylab("") + xlab(""))
        outfile <- sprintf("simulation_methodology_%s_random_initialization_pr_y_%s_sampling_distribution.png", estimator, matrix_index)
        message("saving ", outfile)
        ggsave(outfile, p, width=12.5, height=8)
    }
}

df_melted <- melt(df, id.vars=c("panel_size", "panel_size_label"))
p <- (ggplot(df_melted, aes(x=value)) +
      geom_histogram(binwidth=0.01, fill="white", color="black") +
      facet_wrap(~ variable, scale="free_x"))
p  # Looks correct, compare to params0 -- careful, lots of estimates hit max_iter

## Do min dist estimates vary with initial parameters?
P_hat_range_min_dist <- c(lapply(all_simulations, function(list_of_replications) {
    sapply(list_of_replications, function(replication) {
        matrix_of_P_hat <- sapply(replication$min_dist_params_hat_list, function(x) {
            return(c(x$P_list, recursive=TRUE))
        })
        P_hat_range <- max(apply(matrix_of_P_hat, 1, max) - apply(matrix_of_P_hat, 1, min))
        return(P_hat_range)
    })
}), recursive=TRUE)
sum(P_hat_range_min_dist > 0.01)  # 60 simulations out of length(P_hat_range_min_dist)=1200 where P_hat differs by more than 0.01 across initializations
sum(P_hat_range_min_dist > 0.10)  # 37 simulations out of length(P_hat_range_min_dist)=1200 where P_hat differs by more than 0.10 across initializations
P_hat_range_em <- c(lapply(all_simulations, function(list_of_replications) {
    sapply(list_of_replications, function(replication) {
        matrix_of_P_hat <- sapply(replication$hmm_params_hat_list, function(x) {
            return(c(x$P_list, recursive=TRUE))
        })
        P_hat_range <- max(apply(matrix_of_P_hat, 1, max) - apply(matrix_of_P_hat, 1, min))
        return(P_hat_range)
    })
}), recursive=TRUE)
sum(P_hat_range_em > 0.01)  # 1183 simulations out of length(P_hat_range_em)=1200
sum(P_hat_range_em > 0.10)  # 394 simulations out of length(P_hat_range_em)=1200

## RMSE of each estimator as fn of sample size
dtable <- data.table(df)
setkey(dtable, panel_size)
with(subset(dtable, panel_size <= 5000), t.test(x=(hmm_params_hat_pr_y_11 - params0$pr_y[1, 1])^2,
                                                y=(min_dist_params_hat_pr_y_11 - params0$pr_y[1, 1])^2))  # P-value 0.02
with(subset(dtable, panel_size <= 5000), t.test(x=(hmm_params_hat_P1_11 - params0$P_list[[1]][1, 1])^2,
                                                y=(min_dist_params_hat_P1_11 - params0$P_list[[1]][1, 1])^2))  # P-value 0.11
with(subset(dtable, panel_size <= 5000), t.test(x=(hmm_params_hat_P2_11 - params0$P_list[[2]][1, 1])^2,
                                                y=(min_dist_params_hat_P2_11 - params0$P_list[[2]][1, 1])^2))  # P-value 0.09
with(subset(dtable, panel_size <= 5000), t.test(x=(hmm_params_hat_P3_11 - params0$P_list[[3]][1, 1])^2,
                                                y=(min_dist_params_hat_P3_11 - params0$P_list[[3]][1, 1])^2))  # P-value 0.004
dtable_rmses <- dtable[, list(rmse_hmm_pr_y_11=sqrt(mean((hmm_params_hat_pr_y_11 - params0$pr_y[1, 1])^2)),
                              rmse_min_dist_pr_y_11=sqrt(mean((min_dist_params_hat_pr_y_11 - params0$pr_y[1, 1])^2, na.rm=TRUE)),
                              rmse_hmm_P1_11=sqrt(mean((hmm_params_hat_P1_11 - params0$P_list[[1]][1, 1])^2)),
                              rmse_min_dist_P1_11=sqrt(mean((min_dist_params_hat_P1_11 - params0$P_list[[1]][1, 1])^2, na.rm=TRUE)),
                              rmse_hmm_P2_11=sqrt(mean((hmm_params_hat_P2_11 - params0$P_list[[2]][1, 1])^2)),
                              rmse_min_dist_P2_11=sqrt(mean((min_dist_params_hat_P2_11 - params0$P_list[[2]][1, 1])^2, na.rm=TRUE)),
                              rmse_hmm_P3_11=sqrt(mean((hmm_params_hat_P3_11 - params0$P_list[[3]][1, 1])^2)),
                              rmse_min_dist_P3_11=sqrt(mean((min_dist_params_hat_P3_11 - params0$P_list[[3]][1, 1])^2, na.rm=TRUE))),
                       by="panel_size"]
dtable_melted <- melt(dtable_rmses, id.vars="panel_size")
dtable_melted[, estimator := ifelse(grepl("_hmm_", variable), "EM", "MinDist")]
dtable_melted[, parameter := ifelse(grepl("pr_y_11", variable), "first entry of misclassification matrix",
                                    sprintf("first entry of %s", gsub("rmse_min_dist_|rmse_hmm_|_11", "", variable)))]

p <- (ggplot(dtable_melted, aes(x=panel_size, y=value, color=estimator)) +
      geom_point() + geom_line(size=1.2) +
      scale_color_manual("estimator", values=c("black", "grey")) +
      scale_x_continuous("panel size", breaks=c(100, 1000, 5000, 10000)) +
      scale_y_continuous("root mean squared error") +
      geom_hline(yintercept=0, linetype=2, color="grey") +
      facet_wrap(~ parameter, scales="free_y"))
p
outfile <- "simulation_methodology_random_initialization_rmse_em_and_min_dist.png"
message("saving ", outfile)
ggsave(outfile, p, width=11, height=8)

message("done with simulation, time is ", Sys.time())
