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
