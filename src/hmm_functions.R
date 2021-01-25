get_deforestation_prob_from_P <- function(P) {
    ## Note: this assumes we have 2 hidden states, and state 1 is forest
    return(P[1, 2])
}

get_reforestation_prob_from_P <- function(P) {
    ## Note: this assumes we have 2 hidden states, and state 1 is forest
    return(P[2, 1])
}

get_random_initial_parameters <- function(params0) {

    ## Given a true set of HMM parameters, return random incorrect parameters from which to begin parameter estimation

    initial_parameters <- list(n_components=params0$n_components)

    if ("P_list" %in% names(params0)) {
        initial_parameters$P_list <- lapply(params0$P_list, function(correct_P) {

            ## Probabilities on diagonals of the transition probability matrices
            ## TODO Does min_dist sometimes get stuck at "the wrong" edge of the parameter space, and, if so,
            ## does that happen less frequently if we bump up the minimum value on the diagonals of initial_parameters$P_list?
            random_uniform <- runif(params0$n_components, min=0.60, max=0.98)

            P <- matrix((1 - random_uniform) / (params0$n_components - 1), nrow=nrow(correct_P), ncol=ncol(correct_P))
            diag(P) <- random_uniform

            return(P)
        })
    } else {
        random_uniform <- runif(params0$n_components, min=0.60, max=0.98)

        P <- matrix((1 - random_uniform) / (params0$n_components - 1), nrow=nrow(params0$P), ncol=ncol(params0$P))
        diag(P) <- random_uniform
        initial_parameters$P <- P
    }

    ## Probabilities on the diagonals of the observation probability matrix pr_y
    random_uniform <- runif(params0$n_components, min=0.6, max=0.98)
    initial_parameters$pr_y <- matrix((1 - random_uniform) / (params0$n_components - 1), nrow=nrow(params0$pr_y), ncol=ncol(params0$pr_y))
    diag(initial_parameters$pr_y) <- random_uniform

    initial_parameters$mu <- runif(n=params0$n_components)
    initial_parameters$mu <- initial_parameters$mu / sum(initial_parameters$mu)

    return(initial_parameters)
}

is_diag_dominant <- function(pr_y) {
    return(all(diag(pr_y) > 0.5))
}

get_min_distance_estimates_time_homogeneous <- function(initial_params, M_Y_joint_hat_list, M_Y_joint_hat_inverse_list, M_fixed_y_Y_joint_hat_list, dtable) {

    M_S_joint_initial <- t(initial_params$P * matrix(initial_params$mu, initial_params$n_components, initial_params$n_components))

    x_guess <- c(t(initial_params$pr_y), c(M_S_joint_initial))

    n_components <- initial_params$n_components

    n_equality_constraints <- n_components + 1

    ## Note: the lower bound for the diagonals of pr_y is 0.5 (to make pr_y diagonally dominant);
    ## the lower bound for all other parameters is zero
    ## The index using which(c(diag(initial_params$n_components)) > 0) assumes that the first params$n_components^2 elements
    ## of the x_guess vector (i.e. the argument to objfn_minimum_distance) represent the observation probability matrix pr_y
    lower_bound <- rep(0, length(x_guess))
    lower_bound[which(c(diag(initial_params$n_components)) > 0)] = 0.5

    solnp_result <- solnp(x_guess,
                          fun=objfn_min_dist_time_homogeneous,
                          eqfun=eq_function_min_dist_time_homogeneous,
                          eqB=rep(0, n_equality_constraints),
                          LB=lower_bound,
                          UB=rep(1, length(x_guess)),
                          M_Y_joint_hat_list=M_Y_joint_hat_list,
                          M_Y_joint_hat_inverse_list=M_Y_joint_hat_inverse_list,
                          M_fixed_y_Y_joint_hat_list=M_fixed_y_Y_joint_hat_list,
                          n_components=n_components,
                          control=list(delta=1e-9, tol=1e-13, trace=1, rho=0.1))  # Careful, sensitive to control

    M_Y_given_S_hat_solnp <- matrix(solnp_result$pars[seq(1, n_components^2)], n_components, n_components)  # Transpose of params0$pr_y
    M_S_joint_hat_solnp <- matrix(solnp_result$pars[seq((n_components^2) + 1, (n_components^2)*(2))], n_components, n_components)

    ## Note: we keep track of the objective function values so that we can pick the best MD estimate (lowest objfn_values)
    min_dist_params_hat <- list(pr_y=t(M_Y_given_S_hat_solnp),
                                P=get_transition_probs_from_M_S_joint(M_S_joint_hat_solnp),
                                convergence=solnp_result$convergence,
                                mu=colSums(M_S_joint_hat_solnp),
                                objfn_values=solnp_result$values,
                                x_guess=x_guess)

    return(min_dist_params_hat)

}

get_min_distance_estimates <- function(initial_params, M_Y_joint_hat_list, M_Y_joint_hat_inverse_list, M_fixed_y_Y_joint_hat_list, dtable) {

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

    x_guess <- c(t(initial_params$pr_y), c(M_S_joint_list_initial, recursive=TRUE))
    max_time <- max(dtable$time)

    n_components <- initial_params$n_components

    n_equality_constraints <- n_components + 1 + n_components * (max(dtable$time) - 2)

    ## Note: the lower bound for the diagonals of pr_y is 0.5 (to make pr_y diagonally dominant);
    ## the lower bound for all other parameters is zero
    ## The index using which(c(diag(initial_params$n_components)) > 0) assumes that the first params$n_components^2 elements
    ## of the x_guess vector (i.e. the argument to objfn_minimum_distance) represent the observation probability matrix pr_y
    lower_bound <- rep(0, length(x_guess))
    lower_bound[which(c(diag(initial_params$n_components)) > 0)] = 0.5

    solnp_result <- solnp(x_guess,
                          fun=objfn_minimum_distance,
                          eqfun=eq_function_minimum_distance,
                          eqB=rep(0, n_equality_constraints),
                          LB=lower_bound,
                          UB=rep(1, length(x_guess)),
                          M_Y_joint_hat_list=M_Y_joint_hat_list,
                          M_Y_joint_hat_inverse_list=M_Y_joint_hat_inverse_list,
                          M_fixed_y_Y_joint_hat_list=M_fixed_y_Y_joint_hat_list,
                          max_time=max_time,
                          n_components=n_components,
                          control=list(delta=1e-9, tol=1e-13, trace=1, rho=0.1))  # Careful, sensitive to control

    M_Y_given_S_hat_solnp <- matrix(solnp_result$pars[seq(1, n_components^2)], n_components, n_components)  # Transpose of params0$pr_y
    M_S_joint_list_hat_solnp <- lapply(seq_len(max_time - 1), function(time_index, n_components=initial_params$n_components) {
        return(matrix(solnp_result$pars[seq((n_components^2)*time_index + 1, (n_components^2)*(1 + time_index))], n_components, n_components))

    })
    
    ## Note: we keep track of the objective function values so that we can pick the best MD estimate (lowest objfn_values)
    min_dist_params_hat <- list(pr_y=t(M_Y_given_S_hat_solnp),
                                P_list=lapply(M_S_joint_list_hat_solnp, get_transition_probs_from_M_S_joint),
                                convergence=solnp_result$convergence,
                                mu=colSums(M_S_joint_list_hat_solnp[[1]]),
                                objfn_values=solnp_result$values,
                                x_guess=x_guess)

    return(min_dist_params_hat)

}

get_hmm_and_minimum_distance_estimates_random_initialization <- function(params0, panel, n_random_starts=10) {

    ## Params0 are true HMM parameters used to generate data

    require(data.table)
    require(Rsolnp)

    random_initial_parameters <- replicate(n=n_random_starts, get_random_initial_parameters(params0), simplify=FALSE)

    em_params_hat_list <- lapply(random_initial_parameters, function(initial_params) {
        return(get_expectation_maximization_estimates(panel, initial_params, max_iter=30, epsilon=0.001))
    })
    em_likelihoods <- sapply(em_params_hat_list, function(x) {
        return(max(x$loglik))
    })

    for(idx in seq_along(panel)) {
        panel[[idx]]$point_id <- idx
        panel[[idx]]$time <- seq_along(panel[[idx]]$y)
    }

    dtable <- rbindlist(Map(data.frame, panel))
    setkey(dtable, point_id)

    stopifnot(all(c("point_id", "time", "y") %in% names(dtable)))

    dtable[, y_one_period_ahead := c(tail(y, .N-1), NA), by="point_id"]
    dtable[, y_two_periods_ahead := c(tail(y, .N-2), NA, NA), by="point_id"]

    M_Y_joint_hat_list <- lapply(seq_len(max(dtable$time) - 1), function(fixed_t) {
        with(subset(dtable, time == fixed_t), prop.table(table(y_one_period_ahead, y)))
    })  # Joint distribution of (Y_{t+1}, Y_{t}) and (Y_{t+2}, Y_{t+1})

    ## Compute inverses once and pass them to get_min_distance_estimates / solnp
    M_Y_joint_hat_inverse_list <- lapply(M_Y_joint_hat_list, solve)

    ## TODO Need to handle edge case where any of these matrices are not invertible, which might happen at small sample sizes
    M_fixed_y_Y_joint_hat_list <- lapply(seq_len(params0$n_components), function(fixed_y) {
        lapply(seq_len(max(dtable$time) - 2), function(fixed_t) {
            return(with(subset(dtable, time == fixed_t & y_two_periods_ahead == fixed_y),
                        table(y_one_period_ahead, y)) / sum(dtable$time == fixed_t))
        })
    })

    min_dist_params_hat <- lapply(random_initial_parameters,
                                  get_min_distance_estimates,
                                  M_Y_joint_hat_list=M_Y_joint_hat_list,
                                  M_Y_joint_hat_inverse_list=M_Y_joint_hat_inverse_list,
                                  M_fixed_y_Y_joint_hat_list=M_fixed_y_Y_joint_hat_list,
                                  dtable=dtable)

    min_dist_objfn_values <- sapply(min_dist_params_hat, function(x) {
        return(min(x$objfn_values))
    })

    min_dist_params_hat_best_objfn <- min_dist_params_hat[[which.min(min_dist_objfn_values)]]

    min_dist_pr_y_is_diag_dominant <- sapply(min_dist_params_hat, function(x) {
        return(is_diag_dominant(x$pr_y))
    })

    em_params_hat_best_likelihood <- em_params_hat_list[[which.max(em_likelihoods)]]

    return(list("panel_size"=length(panel),
                "M_Y_joint_hat"=M_Y_joint_hat_list,
                "em_params_hat_list"=em_params_hat_list,
                "em_params_hat_loglikelihoods"=em_likelihoods,
                "initial_parameters_list"=random_initial_parameters,
                "em_params_hat_best_likelihood"=em_params_hat_best_likelihood,
                "min_dist_params_hat"=min_dist_params_hat,
                "min_dist_objfn_values"=min_dist_objfn_values,
                "min_dist_params_hat_best_objfn"=min_dist_params_hat_best_objfn,
                "min_dist_pr_y_is_diag_dominant"=min_dist_pr_y_is_diag_dominant))
}

get_minimum_distance_estimates_random_initialization_time_homogeneous <- function(params0, panel, n_random_starts=10) {

    require(data.table)
    require(Rsolnp)

    random_initial_parameters <- replicate(n=n_random_starts, get_random_initial_parameters(params0), simplify=FALSE)

    ## em_params_hat_list <- lapply(random_initial_parameters, function(initial_params) {
    ##     return(get_expectation_maximization_estimates(panel, initial_params, max_iter=30, epsilon=0.001))
    ## })
    ## em_likelihoods <- sapply(em_params_hat_list, function(x) {
    ##     return(max(x$loglik))
    ## })

    for(idx in seq_along(panel)) {
        panel[[idx]]$point_id <- idx
        panel[[idx]]$time <- seq_along(panel[[idx]]$y)
    }

    dtable <- rbindlist(Map(data.frame, panel))
    setkey(dtable, point_id)

    stopifnot(all(c("point_id", "time", "y") %in% names(dtable)))

    dtable[, y_one_period_ahead := c(tail(y, .N-1), NA), by="point_id"]
    dtable[, y_two_periods_ahead := c(tail(y, .N-2), NA, NA), by="point_id"]

    M_Y_joint_hat_list <- lapply(seq_len(max(dtable$time) - 1), function(fixed_t) {
        with(subset(dtable, time == fixed_t), prop.table(table(y_one_period_ahead, y)))
    })  # Joint distribution of (Y_{t+1}, Y_{t}) and (Y_{t+2}, Y_{t+1})

    ## Compute inverses once and pass them to get_min_distance_estimates / solnp
    M_Y_joint_hat_inverse_list <- lapply(M_Y_joint_hat_list, solve)

    M_fixed_y_Y_joint_hat_list <- lapply(seq_len(params0$n_components), function(fixed_y) {
        lapply(seq_len(max(dtable$time) - 2), function(fixed_t) {
            return(with(subset(dtable, time == fixed_t & y_two_periods_ahead == fixed_y),
                        table(y_one_period_ahead, y)) / sum(dtable$time == fixed_t))
        })
    })

    min_dist_params_hat <- lapply(random_initial_parameters,
                                  get_min_distance_estimates_time_homogeneous,
                                  M_Y_joint_hat_list=M_Y_joint_hat_list,
                                  M_Y_joint_hat_inverse_list=M_Y_joint_hat_inverse_list,
                                  M_fixed_y_Y_joint_hat_list=M_fixed_y_Y_joint_hat_list,
                                  dtable=dtable)

    min_dist_objfn_values <- sapply(min_dist_params_hat, function(x) {
        return(min(x$objfn_values))
    })

    min_dist_params_hat_best_objfn <- min_dist_params_hat[[which.min(min_dist_objfn_values)]]

    min_dist_pr_y_is_diag_dominant <- sapply(min_dist_params_hat, function(x) {
        return(is_diag_dominant(x$pr_y))
    })

    ## em_params_hat_best_likelihood <- em_params_hat_list[[which.max(em_likelihoods)]]

    return(list("panel_size"=length(panel),
                "M_Y_joint_hat"=M_Y_joint_hat_list,
                ## "em_params_hat_list"=em_params_hat_list,
                ## "em_params_hat_loglikelihoods"=em_likelihoods,
                "initial_parameters_list"=random_initial_parameters,
                ## "em_params_hat_best_likelihood"=em_params_hat_best_likelihood,
                "min_dist_params_hat"=min_dist_params_hat,
                "min_dist_objfn_values"=min_dist_objfn_values,
                "min_dist_params_hat_best_objfn"=min_dist_params_hat_best_objfn,
                "min_dist_pr_y_is_diag_dominant"=min_dist_pr_y_is_diag_dominant))
}

get_transition_probs_from_M_S_joint <- function(M_S_joint) {
    return(t(M_S_joint) / rowSums(t(M_S_joint)))
}

objfn_minimum_distance <- function(x, M_Y_joint_hat_inverse_list, M_Y_joint_hat_list, M_fixed_y_Y_joint_hat_list,
                                   max_time, n_components) {
    ## Objective function for minimum distance estimation with time-varying transition probabilities, time-invariant misclassification probabilities
    stopifnot(is.vector(x))
    stopifnot(length(x) == max_time * n_components^2)

    ## TODO Update equation references
    candidate_M_Y_given_S <- matrix(x[seq(1, n_components^2)], n_components, n_components)
    candidate_M_S_joint_list <- lapply(seq_len(max_time - 1), function(time_index) {
        return(matrix(x[seq((n_components^2)*time_index + 1, (n_components^2)*(1 + time_index))], n_components, n_components))
    })
    stopifnot(length(candidate_M_S_joint_list) == length(M_Y_joint_hat_inverse_list))
    candidate_D_list <- lapply(seq_len(n_components), function(fixed_y) {
        lapply(seq_len(max_time - 2), function(fixed_t) {
            candidate_M_S_joint <- candidate_M_S_joint_list[[fixed_t + 1]]
            candidate_P <- t(candidate_M_S_joint / matrix(colSums(candidate_M_S_joint),
                                                          nrow(candidate_M_S_joint),
                                                          ncol(candidate_M_S_joint), byrow=TRUE))  # From t+1 to t+2
            candidate_D <- matrix(0, n_components, n_components)
            diag(candidate_D) <- candidate_P %*% candidate_M_Y_given_S[fixed_y, ]
            return(candidate_D)
        })
    })

    ## TODO This assumes n_components == 3
    if(abs(candidate_M_Y_given_S[1, 1] - candidate_M_Y_given_S[1, 2]) < 0) message('Close to non diag dom')

    stopifnot(length(candidate_D_list) == length(M_fixed_y_Y_joint_hat_list))  # Careful, lists of lists
    stopifnot(length(candidate_D_list[[1]]) == length(M_fixed_y_Y_joint_hat_list[[1]]))  # Careful with fixed_y
    g1_vectors_for_fixed_y <- lapply(seq_along(candidate_D_list), function(fixed_y) {
        sapply(seq_along(candidate_D_list[[fixed_y]]), function(time_index) {
            return(norm(M_fixed_y_Y_joint_hat_list[[fixed_y]][[time_index]] %*%
                        M_Y_joint_hat_inverse_list[[time_index]] %*%
                        candidate_M_Y_given_S -
                        candidate_M_Y_given_S %*%
                        candidate_D_list[[fixed_y]][[time_index]], type="F"))
        })
    })
    g1_vector <- c(g1_vectors_for_fixed_y, recursive=TRUE)
    g2_vector <- sapply(seq_along(candidate_M_S_joint_list), function(time_index) {
        return(norm(candidate_M_Y_given_S %*%
                    candidate_M_S_joint_list[[time_index]] %*% t(candidate_M_Y_given_S) -
                    M_Y_joint_hat_list[[time_index]], type="F"))
    })
    g <- c(g1_vector, g2_vector)

    weights <- NULL
    if(is.null(weights)) {
        weights <- diag(length(g))
    }
    stopifnot(nrow(weights) == length(g) && ncol(weights) == length(g))
    return(as.vector(t(g) %*% weights %*% g))
}

eq_function_minimum_distance <- function(x,
                                         M_Y_joint_hat_inverse_list,
                                         M_Y_joint_hat_list,
                                         M_fixed_y_Y_joint_hat_list,
                                         max_time,
                                         n_components) {
    ## Constraint function for minimum distance estimation (constraint is eq_function(x) = 1 everywhere)
    ## "The main and constraint functions must take the exact same arguments, irrespective of whether they are used"
    candidate_M_Y_given_S <- matrix(x[seq(1, n_components^2)], n_components, n_components)
    candidate_M_S_joint_list <- lapply(seq_len(max_time - 1), function(time_index) {
        return(matrix(x[seq((n_components^2)*time_index + 1, (n_components^2)*(1 + time_index))], n_components, n_components))
    })

    candidate_M_S_rowSums <- sapply(candidate_M_S_joint_list, rowSums)
    candidate_M_S_colSums <- sapply(candidate_M_S_joint_list, colSums)

    differences_in_marginal_distributions <- candidate_M_S_rowSums[, 1:(max_time - 2)] - candidate_M_S_colSums[, 2:(max_time - 1)]

    ## Note: subtract 1 from probabilities so that the contraint function must always equal zero
    return(c(colSums(candidate_M_Y_given_S) - 1.0, sum(candidate_M_S_joint_list[[1]]) - 1.0, differences_in_marginal_distributions))
}

valid_panel_element <- function(panel_element, params) {
    ## Panel element is a list describing a single realization of the HMM
    stopifnot(is.list(panel_element))
    stopifnot("y" %in% names(panel_element))
    stopifnot("pr_y" %in% names(params))
    stopifnot(all(is.na(panel_element$y) | panel_element$y %in% seq_len(ncol(params$pr_y))))  # Discrete in {1, 2, ... , |Y|}
    return(TRUE)
}

valid_parameters <- function(params) {
    stopifnot(is.list(params))
    stopifnot("mu" %in% names(params))  # Vector of probabilities for intial distribution
    stopifnot(length(params$mu) == params$n_components)
    stopifnot(all(params$mu >= 0))
    stopifnot(isTRUE(abs(sum(params$mu)- 1.0)<1e-8))  # TR-- Changed away from Float compairson
    stopifnot(xor("P_list" %in% names(params),  # List of transition matrices for x (one per period)
                  "P" %in% names(params)))  # Time-invariant transition matrix
    stopifnot("pr_y" %in% names(params))  # Observation probabilities conditional on x
    stopifnot("n_components" %in% names(params))
    stopifnot(nrow(params$pr_y) == params$n_components)
    stopifnot(isTRUE(all.equal(rowSums(params$pr_y), rep(1, nrow(params$pr_y)))))  # Float comparison
    return(TRUE)
}

viterbi_path <- function(panel_element, params) {
    ## Viterbi algorithm for HMM with discrete hidden x (returns highest probability path for x)
    ## Written following Ramon van Handel's HMM notes, page 46, algorithm 3.4
    ## https://www.princeton.edu/~rvan/orf557/hmm080728.pdf
    ## Careful, his observation index is in {0, 1, ... , n} while I use {1, 2, ... , t_max}
    stopifnot(valid_panel_element(panel_element, params))
    stopifnot(valid_parameters(params))
    mu <- params$mu
    stopifnot(length(mu) == params$n_components)
    stopifnot(is.list(params$P_list))
    stopifnot(length(params$P_list) == length(panel_element$y) - 1)
    P_list <- params$P_list
    if(is.na(panel_element$y[1])) {
        upsilon <- rep(1, params$n_components)
    } else {
        if("pr_y" %in% names(params)) {
            upsilon <- params$pr_y[, panel_element$y[1]]
        } else {
            upsilon <- params$pr_y_list[[1]][, panel_element$y[1]]
        }
    }
    stopifnot(length(mu) == length(upsilon))  # Same length as state space
    t_max <- length(panel_element$y)
    stopifnot(t_max >= 2)
    v <- matrix(NA, t_max, params$n_components)  # Maximized log likelihoods
    v[1, ] <- log(mu) + log(upsilon)
    b <- matrix(NA, t_max, params$n_components)
    for(k in seq(2, t_max)) {
        P <- P_list[[k - 1]]
        if(is.na(panel_element$y[k])) {
            upsilon <- rep(1, params$n_components)
        } else {
            if("pr_y" %in% names(params)) {
                upsilon <- params$pr_y[, panel_element$y[k]]
            } else {
                upsilon <- params$pr_y_list[[k]][, panel_element$y[k]]
            }
        }
        stopifnot(length(upsilon) == ncol(v))
        for(i in seq_len(params$n_components)) {
            b[k, i] <- which.max(v[k-1, ] + log(P[, i]))
            v[k, i] <- v[k-1, b[k, i]] + log(P[b[k, i], i]) + log(upsilon[i])
        }
    }
    stopifnot(all(!is.na(v)))  # Entire v matrix should be populated
    most_likely_path <- rep(NA, t_max)
    most_likely_path[t_max] <- which.max(v[t_max, ])
    for(k in seq(1, t_max - 1)) {
        most_likely_path[t_max - k] <- b[t_max - k + 1, most_likely_path[t_max - k + 1]]
    }
    stopifnot(is.vector(most_likely_path))
    stopifnot(length(most_likely_path) == t_max)
    stopifnot(all(!is.na(most_likely_path)))
    stopifnot(all(most_likely_path %in% seq_len(params$n_components)))
    return(most_likely_path)
}

apply_viterbi_path_in_parallel <- function(panel, params_hat, max_cores=30) {

    ## Apply viterbi to every element in panel
    num_cores <- min(detectCores(), max_cores)
    cluster <- makeCluster(num_cores)  # Call stopCluster when done

    vars_to_export <- c("viterbi_path", "valid_panel_element", "valid_parameters")
    clusterExport(cl=cluster, varlist=vars_to_export, envir=.GlobalEnv)

    list_of_viterbi_paths <- parLapply(cluster, panel, viterbi_path, params=params_hat)
    stopCluster(cluster)

    return(list_of_viterbi_paths)
}

baum_welch <- function(panel_element, params) {
    ## Baum-Welch algorithm for HMM with discrete hidden x, discrete observations y (NAs allowed)
    ## Written following Ramon van Handel's HMM notes, page 40, algorithm 3.2
    ## https://www.princeton.edu/~rvan/orf557/hmm080728.pdf
    ## Careful, his observation index is in {0, 1, ... , n} while I use {1, 2, ... , y_length}
    stopifnot(valid_panel_element(panel_element, params))
    stopifnot(valid_parameters(params))
    y_length <- length(panel_element$y)
    stopifnot(y_length > 1)
    c <- vector("numeric", y_length)
    if(is.na(panel_element$y[1])) {
        upsilon <- rep(1, params$n_components)
    } else {
        if("pr_y" %in% names(params)) {
            upsilon <- params$pr_y[, panel_element$y[1]]
        } else {
            upsilon <- params$pr_y_list[[1]][, panel_element$y[1]]
        }
    }
    stopifnot(length(upsilon) == params$n_components)
    mu <- params$mu
    c[1] <- sum(upsilon * mu)
    ## Matrix pi_contemporaneous gives probabilities over x_k conditional on {y_1, y_2, ... , y_k}
    ## Notation in van Handel's HMM notes is pi_k, whereas pi_{k|n} conditions on full history of y
    pi_contemporaneous <- matrix(NA, params$n_components, y_length)
    pi_contemporaneous[, 1] <- upsilon * mu / c[1]
    P_list <- params$P_list
    P_transpose_list <- lapply(P_list, t)
    upsilon_list <- list()
    upsilon_list[[1]] <- upsilon
    for(k in seq(2, y_length)) {
        ## Forward loop
        if(is.na(panel_element$y[k])) {
            upsilon <- rep(1, params$n_components)
        } else {
            if("pr_y" %in% names(params)) {
                upsilon <- params$pr_y[, panel_element$y[k]]
            } else {
                upsilon <- params$pr_y_list[[k]][, panel_element$y[k]]
            }
        }
        upsilon_list[[k]] <- upsilon  # Cache for backward loop
        pi_tilde <- upsilon * P_transpose_list[[k-1]] %*% pi_contemporaneous[, k-1]
        c[k] <- sum(pi_tilde)
        pi_contemporaneous[, k] <- pi_tilde  / c[k]
    }
    beta <- matrix(NA, params$n_components, y_length)
    beta[, y_length] <- 1 / c[y_length]
    ## Matrix pi gives probabilities over x conditional on full history of y
    ## Notation in van Handel's HMM notes is pi_{k|n}, as opposed to pi_k
    pi <- matrix(NA, params$n_components, y_length)
    pi[, y_length] <- pi_contemporaneous[, y_length]
    pi_transition_list <- list()  # List of posterior probabilities over hidden x transitions
    for(k in seq(1, y_length - 1)) {
        ## Backward loop
        upsilon <- diag(upsilon_list[[y_length - k + 1]], params$n_components, params$n_components)
        pi_matrix <- diag(pi_contemporaneous[, y_length - k],
                          params$n_components, params$n_components)
        beta_matrix <- diag(beta[, y_length - k + 1], params$n_components, params$n_components)
        P_matrix <- P_list[[y_length - k]]
        beta[, y_length - k] <- P_matrix %*% upsilon %*% beta[, y_length - k + 1] / c[y_length - k]
        pi_transition_list[[y_length - k]] <- pi_matrix %*% P_matrix %*% upsilon %*% beta_matrix
        stopifnot(isTRUE(all.equal(sum(pi_transition_list[[y_length - k]]), 1.0)))
        pi[, y_length - k] <- rowSums(pi_transition_list[[y_length - k]])
    }
    loglik <- sum(log(c))
    return(list(loglik=loglik, pi=pi, pi_transition_list=pi_transition_list))
}

get_expectation_maximization_estimates <- function(panel, params, max_iter, epsilon=0.001) {
    ## EM for panel of independent HMM realizations; stop at max_iter or distance < epsilon
    ## Written following Ramon van Handel's HMM notes page 87, algorithm 6.1, modified for panel data
    ## https://www.princeton.edu/~rvan/orf557/hmm080728.pdf
    message("starting em ", Sys.time())
    stopifnot(valid_parameters(params))
    observation_lengths <- vapply(panel, function(panel_element) {
        length(panel_element$y)
    }, FUN.VALUE=1)  # Observation lengths should be the same in each panel element
    stopifnot(all(observation_lengths > 1) && length(unique(observation_lengths)) == 1)
    observation_length <- observation_lengths[1]
    iteration <- 1
    while(iteration <= max_iter) {
        baum_welch_list <- lapply(panel, baum_welch, params=params)
        panel_with_baum_welch <- Map(c, panel, baum_welch_list)
        ## Update list of matrices containing transition probabilities (one per timeperiod)
        updated_P_list <- lapply(seq(1, observation_length - 1), function(time) {
            numerators <- lapply(baum_welch_list, function(baum_welch_element) {
                baum_welch_element$pi_transition_list[[time]]
            })
            denominators <- lapply(baum_welch_list, function(baum_welch_element) {
                matrix(baum_welch_element$pi[, time], params$n_components, params$n_components)
            })
            numerator <- Reduce("+", numerators)
            denominator <- Reduce("+", denominators)
            stopifnot(all(rowSums(denominator) > 0))  # Can fail if a state has zero probability
            stopifnot(all(denominator > 0))
                P <- numerator / denominator
            stopifnot(isTRUE(all.equal(rowSums(P), rep(1, params$n_components))))
            return(P)
        })
        ## Update initial probabilities over hidden state
        mus <- lapply(baum_welch_list, function(baum_welch_element) {
            baum_welch_element$pi[, 1]  # Initial distribution
        })
        updated_mu <- Reduce("+", mus) / length(panel)
        stopifnot(isTRUE(all.equal(sum(updated_mu), 1)))
        ## Update observation probabilities; careful, observation vector can contain NAs
        pr_y_weights <- lapply(panel_with_baum_welch, function(panel_element, params) {
            y_non_NA <- matrix(1 * !is.na(panel_element$y),
                               nrow(params$pr_y), length(panel_element$y), byrow=TRUE)
            rowSums(panel_element$pi * y_non_NA)  # Sum over time
        }, params)
        pr_y_matrices <- lapply(panel_with_baum_welch, function(panel_element, params) {
            column_list <- lapply(seq_len(ncol(params$pr_y)), function(curr_y) {
                y_indicators <- matrix(!is.na(panel_element$y) &
                                       panel_element$y == curr_y,
                                       nrow(params$pr_y),
                                       length(panel_element$y), byrow=TRUE)
                return(rowSums(panel_element$pi * y_indicators))  # Sum over time
            })
            return(do.call(cbind, column_list))
        }, params)
        updated_pr_y <- (Reduce("+", pr_y_matrices) / Reduce("+", pr_y_weights))
        stopifnot(isTRUE(all.equal(rowSums(updated_pr_y), rep(1, nrow(updated_pr_y)))))
        ## Take sup norm between previous parameters and updated values
        if("P_list" %in% names(params)) {
            P_distances <- vapply(seq_along(params$P_list), function(i) {
                max(abs(params$P_list[[i]] - updated_P_list[[i]]))
            }, FUN.VALUE=1)
        } else {
            P_distances <- vapply(seq_along(params$P_coef_list), function(i) {
                max(abs(c(params$P_coef_list[[i]], recursive=TRUE) -
                        c(updated_P_coef_list[[i]], recursive=TRUE)))
            }, FUN.VALUE=1)
        }
        if("mu" %in% names(params)) {
            mu_distance <- max(abs(params$mu - updated_mu))
        } else if("mu_coefs" %in% names(params)) {
            mu_distance <- max(abs(c(params$mu_coefs, recursive=TRUE) -
                                   c(updated_mu_coefs, recursive=TRUE)))
        }
        if("pr_y" %in% names(params)) {
            pr_y_distance <- max(abs(params$pr_y - updated_pr_y))
        } else {
            pr_y_distance <- max(abs(c(params$pr_y_list, recursive=TRUE) - c(updated_pr_y_list, recursive=TRUE)))
        }
        distance <- max(mu_distance, P_distances, pr_y_distance)
        loglik <- sum(vapply(baum_welch_list, function(x) x$loglik, FUN.VALUE=1))
        message("iteration ", iteration,
                " distance ", round(distance, 4), " loglik ", round(loglik, 4))
        if("distance" %in% names(params)) {
            params$distance <- c(params$distance, distance)
        } else {
            params$distance <- distance  # Maximum distance over all parameters
        }
        if("mu_distance" %in% names(params)) {
            params$mu_distance <- c(params$mu_distance, mu_distance)
        } else {
            params$mu_distance <- mu_distance  # Distance for initial distribution over x
        }
        if("pr_y_distance" %in% names(params)) {
            params$pr_y_distance <- c(params$pr_y_distance, pr_y_distance)
        } else {
            params$pr_y_distance <- pr_y_distance
        }
        if("P_distance" %in% names(params)) {
            params$P_distance <- c(params$P_distance, max(P_distances))
        } else {
            params$P_distance <- max(P_distances)
        }
        if("n_em_iterations" %in% names(params)) {
            params$n_em_iterations <- params$n_em_iterations + 1  # Could be useful to track
        } else {
            params$n_em_iterations <- 1
        }
        if("mu" %in% names(params)) {
            params$mu <- updated_mu
        } else if("mu_coefs" %in% names(params)) {
            params$mu_coefs <- updated_mu_coefs
        }
        if("P_list" %in% names(params)) {
            params$P_list <- updated_P_list
        } else {
            params$P_coef_list <- updated_P_coef_list
        }
        if("pr_y" %in% names(params)) {
            params$pr_y <- updated_pr_y
        } else {
            params$pr_y_list <- updated_pr_y_list
        }
        if("loglik" %in% names(params)) {
            params$loglik <- c(params$loglik, loglik)  # Save history of log likelihoods
        } else {
            params$loglik <- loglik
        }
        params$panel_size <- length(panel)
        if(distance < epsilon) break
        iteration <- iteration + 1
        rm(panel_with_baum_welch)
        rm(baum_welch_list)
        gc()
    }
    message("done running em ", Sys.time())
    params$time_finished_em <- Sys.time()
    return(params)
}

simulate_discrete_markov <- function(params) {
    ## Simulate distrete markov chain with transitions P_list, initial distribution mu
    stopifnot(valid_parameters(params))
    stopifnot(all(c("mu",
                    "P_list") %in% names(params)))  # Does not accept mu_coefs or P_coef_list
    x_length <- length(params$P_list) + 1
    stopifnot(x_length > 1)
    stopifnot(all(vapply(params$P_list, function(P) {
        isTRUE(all.equal(rowSums(P), rep(1, nrow(P))))
    }, FUN.VALUE=TRUE)))
    params$P_list_dims <- vapply(params$P_list, dim, FUN.VALUE=c(0, 1))
    stopifnot(length(unique(as.vector(params$P_list_dims))) == 1)  # Same nrow, ncol in all P
    state_space <- seq_len(params$n_components)
    x <- vector("numeric", x_length)
    x[1] <- sample(state_space, 1, prob=params$mu)
    for (t in seq(2, x_length)) {
        x[t] <- sample(state_space, 1, prob=params$P_list[[t - 1]][x[t - 1], ])
    }
    stopifnot(all(x %in% state_space))
    return(x)
}

simulate_hmm <- function(params) {
    stopifnot(valid_parameters(params))
    ## TODO This is called S in the paper, not X.  Call it "state" for clarity?
    x <- simulate_discrete_markov(params)
    if("pr_y" %in% names(params)) {
        y <- vapply(x, function(x) {
            sample(seq_len(ncol(params$pr_y)), size=1, prob=params$pr_y[x, ])
        }, FUN.VALUE=1)
    } else {
        y <- vector("numeric", length=length(x))
        for(time in seq_along(y)) {
            y[time] <- sample(seq_len(ncol(params$pr_y_list[[time]])), size=1, prob=params$pr_y_list[[time]][x[time], ])
        }
    }
    return(list(x=x, y=y))
}

rows_sum_to_one <- function(probability_matrix) {
    ## Useful for checking transition and observation probabilities
    return(isTRUE(all.equal(rowSums(probability_matrix), rep(1, nrow(probability_matrix)))))
}

get_hmm_panel_from_points <- function(points_dt, discrete_y_varname, max_panel_size=NULL) {
    ## Given a points_dt data.table, return a panel in format expected by EM estimation function
    ## If max_panel_size is non-NULL and max_panel_size < length(unique(points_dt$point_id)), return a random sample
    stopifnot(is.data.table(points_dt))
    stopifnot("point_id" %in% names(points_dt))
    stopifnot(discrete_y_varname %in% names(points_dt))
    if(!is.null(max_panel_size) && max_panel_size < length(unique(points_dt$point_id))) {
        message("generating panel from points_dt (taking sample of size ", max_panel_size, ")")
        sample_point_id <- sample(unique(points_dt$point_id), size=max_panel_size)
    } else {
        message("generating panel from points_dt (using full sample, keeping order unchanged)")
        sample_point_id <- unique(points_dt$point_id)
    }
    panel <- lapply(sample_point_id, function(curr_point_id) {
        curr_rows <- points_dt[J(curr_point_id)]
        discrete_y <- as.integer(curr_rows[[discrete_y_varname]])  # From factor to integer
        panel_element <- list(point_id=curr_point_id, y=discrete_y)
        if("validation_landuse" %in% names(points_dt)) {
            panel_element$validation_landuse <- points_dt[J(curr_point_id)]$validation_landuse
        }
        if("validation_landuse_coarse" %in% names(points_dt)) {
            panel_element$validation_landuse_coarse <- points_dt[J(curr_point_id)]$validation_landuse_coarse
        }
        return(panel_element)
    })
    return(panel)
}

baum_welch_time_homogeneous <- function(panel_element, params) {
    ## Baum-Welch algorithm for HMM with discrete hidden x, discrete observations y (NAs allowed)
    ## Written following Ramon van Handel's HMM notes, page 40, algorithm 3.2
    ## https://www.princeton.edu/~rvan/orf557/hmm080728.pdf
    ## Careful, his observation index is in {0, 1, ... , n} while I use {1, 2, ... , y_length}
    stopifnot(valid_panel_element(panel_element, params))
    stopifnot(valid_parameters_time_homogeneous(params))
    y_length <- length(panel_element$y)
    stopifnot(y_length > 1)
    c <- vector("numeric", y_length)
    if(is.na(panel_element$y[1])) {
        upsilon <- rep(1, params$n_components)
    } else {
        upsilon <- params$pr_y[, panel_element$y[1]]  # Time-invariant pr_y
    }
    stopifnot(length(upsilon) == params$n_components)
    mu <- params$mu
    c[1] <- sum(upsilon * mu)
    ## Matrix pi_contemporaneous gives probabilities over x_k conditional on {y_1, y_2, ... , y_k}
    ## Notation in van Handel's HMM notes is pi_k, whereas pi_{k|n} conditions on full history of y
    pi_contemporaneous <- matrix(NA, params$n_components, y_length)
    pi_contemporaneous[, 1] <- upsilon * mu / c[1]
    upsilon_list <- list()
    upsilon_list[[1]] <- upsilon
    for(k in seq(2, y_length)) {
        ## Forward loop
        if(is.na(panel_element$y[k])) {
            upsilon <- rep(1, params$n_components)
        } else {
            upsilon <- params$pr_y[, panel_element$y[k]]  # Time-invariant pr_y
        }
        upsilon_list[[k]] <- upsilon  # Cache for backward loop
        pi_tilde <- upsilon * t(params$P) %*% pi_contemporaneous[, k-1]
        c[k] <- sum(pi_tilde)
        pi_contemporaneous[, k] <- pi_tilde  / c[k]
    }
    beta <- matrix(NA, params$n_components, y_length)
    beta[, y_length] <- 1 / c[y_length]
    ## Matrix pi gives probabilities over x conditional on full history of y
    ## Notation in van Handel's HMM notes is pi_{k|n}, as opposed to pi_k
    pi <- matrix(NA, params$n_components, y_length)
    pi[, y_length] <- pi_contemporaneous[, y_length]
    pi_transition_list <- list()  # List of posterior probabilities over hidden x transitions
    for(k in seq(1, y_length - 1)) {
        ## Backward loop
        upsilon <- diag(upsilon_list[[y_length - k + 1]], params$n_components, params$n_components)
        pi_matrix <- diag(pi_contemporaneous[, y_length - k],
                          params$n_components, params$n_components)
        beta_matrix <- diag(beta[, y_length - k + 1], params$n_components, params$n_components)
        beta[, y_length - k] <- params$P %*% upsilon %*% beta[, y_length - k + 1] / c[y_length - k]
        pi_transition_list[[y_length - k]] <- pi_matrix %*% params$P %*% upsilon %*% beta_matrix
        stopifnot(isTRUE(all.equal(sum(pi_transition_list[[y_length - k]]), 1.0)))
        pi[, y_length - k] <- rowSums(pi_transition_list[[y_length - k]])
    }
    loglik <- sum(log(c))
    return(list(loglik=loglik, pi=pi, pi_transition_list=pi_transition_list))
}

em_parameter_estimates_time_homogeneous <- function(panel, params, max_iter, epsilon=0.001) {
    ## EM for panel of independent HMM realizations; stop at max_iter or distance < epsilon
    ## Written following Ramon van Handel's HMM notes page 87, algorithm 6.1, modified for panel data
    ## https://www.princeton.edu/~rvan/orf557/hmm080728.pdf
    message("starting em ", Sys.time())
    stopifnot(valid_parameters_time_homogeneous(params))
    observation_lengths <- vapply(panel, function(panel_element) {
        length(panel_element$y)
    }, FUN.VALUE=1)  # Observation lengths should be the same in each panel element
    stopifnot(all(observation_lengths > 1) && length(unique(observation_lengths)) == 1)
    observation_length <- observation_lengths[1]
    iteration <- 1
    while(iteration <= max_iter) {
        baum_welch_list <- lapply(panel, baum_welch_time_homogeneous, params=params)
        panel_with_baum_welch <- Map(c, panel, baum_welch_list)
        ## Update transition probabilities
        numerators <- unlist(lapply(baum_welch_list, function(baum_welch_element) {
            return(baum_welch_element$pi_transition_list)
        }), recursive=FALSE)
        stopifnot(is.list(numerators) && is.matrix(numerators[[1]]))
        denominators <- unlist(lapply(baum_welch_list, function(baum_welch_element) {
            lapply(seq_len(ncol(baum_welch_element$pi) - 1), function(time) {
                matrix(baum_welch_element$pi[, time], params$n_components, params$n_components)
            })
        }), recursive=FALSE)
        stopifnot(is.list(denominators) && is.matrix(denominators[[1]]))
        stopifnot(length(numerators) == length(denominators))
        numerator <- Reduce("+", numerators)
        denominator <- Reduce("+", denominators)
        stopifnot(all(rowSums(denominator) > 0))  # Can fail if a state has zero probability
        stopifnot(all(denominator > 0))
        updated_P <- numerator / denominator
        stopifnot(isTRUE(all.equal(rowSums(updated_P), rep(1, params$n_components))))
        ## Update initial probabilities over hidden state
        mus <- lapply(baum_welch_list, function(baum_welch_element) {
            baum_welch_element$pi[, 1]  # Initial distribution
        })
        updated_mu <- Reduce("+", mus) / length(panel)
        stopifnot(isTRUE(all.equal(sum(updated_mu), 1)))
        ## Update time-invariant observation probabilities; careful, observation vector can contain NAs
        pr_y_weights <- lapply(panel_with_baum_welch, function(panel_element, params) {
            y_non_NA <- matrix(1 * !is.na(panel_element$y),
                               nrow(params$pr_y), length(panel_element$y), byrow=TRUE)
            rowSums(panel_element$pi * y_non_NA)  # Sum over time
        }, params)
        pr_y_matrices <- lapply(panel_with_baum_welch, function(panel_element, params) {
            column_list <- lapply(seq_len(ncol(params$pr_y)), function(curr_y) {
                    y_indicators <- matrix(!is.na(panel_element$y) &
                                           panel_element$y == curr_y,
                                           nrow(params$pr_y),
                                           length(panel_element$y), byrow=TRUE)
                    return(rowSums(panel_element$pi * y_indicators))  # Sum over time
            })
            return(do.call(cbind, column_list))
        }, params)
        updated_pr_y <- (Reduce("+", pr_y_matrices) / Reduce("+", pr_y_weights))
        stopifnot(isTRUE(all.equal(rowSums(updated_pr_y), rep(1, nrow(updated_pr_y)))))
        ## Take sup norm between previous parameters and updated values
        P_distance <- max(abs(params$P - updated_P))
        mu_distance <- max(abs(params$mu - updated_mu))
        pr_y_distance <- max(abs(params$pr_y - updated_pr_y))
        distance <- max(mu_distance, P_distance, pr_y_distance)
        loglik <- sum(vapply(baum_welch_list, function(x) x$loglik, FUN.VALUE=1))
        message("iteration ", iteration, " distance ", round(distance, 4), " loglik ", round(loglik, 4))
        if("distance" %in% names(params)) {
            params$distance <- c(params$distance, distance)
        } else {
            params$distance <- distance  # Maximum distance over all parameters
        }
        if("mu_distance" %in% names(params)) {
            params$mu_distance <- c(params$mu_distance, mu_distance)
        } else {
            params$mu_distance <- mu_distance  # Distance for initial distribution over x
        }
        if("pr_y_distance" %in% names(params)) {
            params$pr_y_distance <- c(params$pr_y_distance, pr_y_distance)
        } else {
            params$pr_y_distance <- pr_y_distance
        }
        if("P_distance" %in% names(params)) {
            params$P_distance <- c(params$P_distance, max(P_distance))
        } else {
            params$P_distance <- max(P_distance)
        }
        if("n_em_iterations" %in% names(params)) {
            params$n_em_iterations <- params$n_em_iterations + 1  # Could be useful to track
        } else {
            params$n_em_iterations <- 1
        }
        params$mu <- updated_mu
        params$P <- updated_P
        params$pr_y <- updated_pr_y
        if("loglik" %in% names(params)) {
            params$loglik <- c(params$loglik, loglik)  # Save history of log likelihoods
        } else {
            params$loglik <- loglik
        }
        params$panel_size <- length(panel)
        if(distance < epsilon) break
        iteration <- iteration + 1
        rm(panel_with_baum_welch)
        rm(baum_welch_list)
        gc()
    }
    message("done running em ", Sys.time())
    params$time_finished_em <- Sys.time()
    return(params)
}

objfn_min_dist_time_homogeneous <- function(x, M_Y_joint_hat_inverse_list, M_Y_joint_hat_list, M_fixed_y_Y_joint_hat_list,
                                            n_components, W_matrix=NULL) {
    ## Objective function for minimum distance estimation
    stopifnot(is.vector(x))
    stopifnot(length(x) == 2 * n_components^2)  # Misclassification probabilities and joint distribution

    candidate_M_Y_given_S <- matrix(x[seq(1, n_components^2)], n_components, n_components)
    candidate_M_S_joint <- matrix(x[seq((n_components^2) + 1, (n_components^2)*(2))], n_components, n_components)

    candidate_D_list <- lapply(seq_len(n_components), function(fixed_y) {
        lapply(seq_len(length(unique(dtable$year)) - 2), function(fixed_t) {
            candidate_P <- candidate_M_S_joint / matrix(colSums(candidate_M_S_joint),
                                                        nrow(candidate_M_S_joint),
                                                        ncol(candidate_M_S_joint), byrow=TRUE)  # From t+1 to t+2
            candidate_D <- matrix(0, n_components, n_components)
            diag(candidate_D) <- candidate_P %*% candidate_M_Y_given_S[fixed_y, ]
            return(candidate_D)
        })
    })

    g1_vectors_for_fixed_y <- lapply(seq_along(candidate_D_list), function(fixed_y) {
        sapply(seq_along(candidate_D_list), function(time_index) {
            return(norm(M_fixed_y_Y_joint_hat_list[[fixed_y]][[time_index]] %*% M_Y_joint_hat_inverse_list[[time_index]] %*%
                        candidate_M_Y_given_S -
                        candidate_M_Y_given_S %*% candidate_D_list[[fixed_y]][[time_index]], type="F"))
        })
    })
    g1_vector <- c(g1_vectors_for_fixed_y, recursive=TRUE)

    g2_vector <- sapply(seq_along(M_Y_joint_hat_list), function(time_index) {
        return(norm(candidate_M_Y_given_S %*% candidate_M_S_joint %*% t(candidate_M_Y_given_S) -
                    M_Y_joint_hat_list[[time_index]], type="F"))
    })

    g <- c(g1_vector, g2_vector)
    if(is.null(W_matrix)) {
        W_matrix <- diag(length(g))
    }
    stopifnot(nrow(W_matrix) == length(g) && ncol(W_matrix) == length(g))

    return(as.vector(t(g) %*% W_matrix %*% g))
}

eq_function_min_dist_time_homogeneous <- function(x, M_Y_joint_hat_inverse_list, M_Y_joint_hat_list, M_fixed_y_Y_joint_hat_list,
                                                  n_components, W_matrix=NULL) {
    ## Constraint function for minimum distance estimation (constraint is eq_function(x) = 1 everywhere)
    ## "The main and constraint functions must take the exact same arguments, irrespective of whether they are used"

    candidate_M_Y_given_S <- matrix(x[seq(1, n_components^2)], n_components, n_components)
    candidate_M_S_joint <- matrix(x[seq((n_components^2) + 1, (n_components^2)*(2))], n_components, n_components)

    ## Note: subtract 1 from probabilities so that the contraint function must always equal zero
    return(c(colSums(candidate_M_Y_given_S) - 1.0, sum(candidate_M_S_joint) - 1.0))
}

valid_parameters_time_homogeneous <- function(params) {
    stopifnot(is.list(params))
    stopifnot("n_components" %in% names(params))
    stopifnot("mu" %in% names(params))  # Vector of probabilities for intial distribution
    stopifnot(length(params$mu) == params$n_components)
    stopifnot(all(params$mu >= 0))
    stopifnot(isTRUE(all.equal(sum(params$mu), 1.0)))  # Careful with float comparisons
    stopifnot("P" %in% names(params))  # Transition matrix for hidden state
    stopifnot("pr_y" %in% names(params))  # Observation probabilities conditional on x
    stopifnot(nrow(params$pr_y) == params$n_components)
    stopifnot(isTRUE(all.equal(rowSums(params$pr_y), rep(1, nrow(params$pr_y)))))  # Float comparison
    return(TRUE)
}

viterbi_path_time_homogeneous <- function(panel_element, params) {
    ## Viterbi algorithm for HMM with discrete hidden x (returns highest probability path for x)
    ## Written following Ramon van Handel's HMM notes, page 46, algorithm 3.4
    ## https://www.princeton.edu/~rvan/orf557/hmm080728.pdf
    ## Careful, his observation index is in {0, 1, ... , n} while I use {1, 2, ... , t_max}
    stopifnot(valid_panel_element(panel_element, params))
    stopifnot(valid_parameters_time_homogeneous(params))
    mu <- params$mu
    stopifnot(length(mu) == params$n_components)
    if(is.na(panel_element$y[1])) {
        upsilon <- rep(1, params$n_components)
    } else {
        if("pr_y" %in% names(params)) {
            upsilon <- params$pr_y[, panel_element$y[1]]
        } else {
            upsilon <- params$pr_y_list[[1]][, panel_element$y[1]]
        }
    }
    stopifnot(length(mu) == length(upsilon))  # Same length as state space
    t_max <- length(panel_element$y)
    stopifnot(t_max >= 2)
    v <- matrix(NA, t_max, params$n_components)  # Maximized log likelihoods
    v[1, ] <- log(mu) + log(upsilon)
    b <- matrix(NA, t_max, params$n_components)
    for(k in seq(2, t_max)) {
        P <- params$P
        if(is.na(panel_element$y[k])) {
            upsilon <- rep(1, params$n_components)
        } else {
            if("pr_y" %in% names(params)) {
                upsilon <- params$pr_y[, panel_element$y[k]]
            } else {
                upsilon <- params$pr_y_list[[k]][, panel_element$y[k]]
            }
        }
        stopifnot(length(upsilon) == ncol(v))
        for(i in seq_len(params$n_components)) {
            b[k, i] <- which.max(v[k-1, ] + log(P[, i]))
            v[k, i] <- v[k-1, b[k, i]] + log(P[b[k, i], i]) + log(upsilon[i])
        }
    }
    stopifnot(all(!is.na(v)))  # Entire v matrix should be populated
    most_likely_path <- rep(NA, t_max)
    most_likely_path[t_max] <- which.max(v[t_max, ])
    for(k in seq(1, t_max - 1)) {
        most_likely_path[t_max - k] <- b[t_max - k + 1, most_likely_path[t_max - k + 1]]
    }
    stopifnot(is.vector(most_likely_path))
    stopifnot(length(most_likely_path) == t_max)
    stopifnot(all(!is.na(most_likely_path)))
    stopifnot(all(most_likely_path %in% seq_len(params$n_components)))
    return(most_likely_path)
}
