rm(list=ls())
library(Rsolnp)
library(data.table)
library(ggplot2)
library(grid)
library(parallel)
library(plyr)
library(stringr)

## Run MC simulations using regression

source("hmm_functions.R")
source("hmm_parameters.R")

set.seed(321321)

## We will run regressions with cell-level data
## Cells could be regions in a country, for example
nCells <- 100

## These are predictors / Xs / covariates that vary at the cell level. These X's are binary -- simulating a policy implemented in 20% of cells/regions
xVec <- runif(nCells) < .2

## The true deforestation rate (in period 3) varies with X, at the cell level
alpha <- .15
beta <- -.10
trueDeforestPr <- alpha + xVec * beta

nPixelPerCell <- 1000

nSim <- 12  # TODO Turn back up to 100

params0 <- get_params0()

## Initially, params0 has 3 transition probability matrices (i.e. 4 periods)
## Let's add more periods (so that we simulate longer timeseries at each pixel)
n_extra_periods <- 3
for(i in seq_len(n_extra_periods)) {
    params0$P_list[[length(params0$P_list) + 1]] <- params0$P_list[[1]]
}

tMax <- length(params0$P_list) + 1  # Number of time periods

## Deforestation only varies (with respect to xVec) in the 3rd period
## These parameters are held fixed across simulations
deforestation_time_index <- 3  # This is the period in which we are interested in measuring how deforestation varies wrt xVec
paramList <- list()
for (i in 1:length(xVec)) {
    paramList[[i]] <- params0
    paramList[[i]]$P_list[[deforestation_time_index]][1, 2] <- trueDeforestPr[i]
    paramList[[i]]$P_list[[deforestation_time_index]][1, 1] <- 1 - trueDeforestPr[i]
}

simulate_one_cell <- function(params) {
    panel <- replicate(nPixelPerCell, simulate_hmm(params), simplify=FALSE)
    for(idx in seq_along(panel)) {
        panel[[idx]]$point_id <- idx
        panel[[idx]]$time <- seq_along(panel[[idx]]$y)
    }
    return(panel)
}

get_estimated_deforestation_rate_viterbi <- function(panel) {
    dtable <- rbindlist(Map(data.frame, panel))
    setkey(dtable, point_id)
    dtable[, y_viterbi_one_period_ahead := c(tail(y_viterbi, .N-1), NA), by="point_id"]

    M_Y_joint_hat_list <- lapply(seq_len(max(dtable$time) - 1), function(fixed_t) {
        with(subset(dtable, time == fixed_t), prop.table(table(y_viterbi_one_period_ahead, y)))
    })

    return(get_transition_probs_from_M_S_joint(M_Y_joint_hat_list[[deforestation_time_index]])[1, 2])
}

get_estimated_deforestation_rate_observations <- function(panel) {
    dtable <- rbindlist(Map(data.frame, panel))
    setkey(dtable, point_id)
    dtable[, y_one_period_ahead := c(tail(y, .N-1), NA), by="point_id"]

    M_Y_joint_hat_list <- lapply(seq_len(max(dtable$time) - 1), function(fixed_t) {
        with(subset(dtable, time == fixed_t), prop.table(table(y_one_period_ahead, y)))
    })

    return(get_transition_probs_from_M_S_joint(M_Y_joint_hat_list[[deforestation_time_index]])[1, 2])
}

get_estimated_deforestation_rate_ground_truth <- function(panel) {
    dtable <- rbindlist(Map(data.frame, panel))
    setkey(dtable, point_id)

    ## Note: we use S for the hidden state in the paper, but it's called X in the simulate code
    ## TODO Would be less confusing to use S in the code, can clean this up later
    dtable[, x_one_period_ahead := c(tail(x, .N-1), NA), by="point_id"]

    M_S_joint_hat_list <- lapply(seq_len(max(dtable$time) - 1), function(fixed_t) {
        with(subset(dtable, time == fixed_t), prop.table(table(x_one_period_ahead, x)))
    })

    return(get_transition_probs_from_M_S_joint(M_S_joint_hat_list[[deforestation_time_index]])[1, 2])
}

regResList <- list()
for (s in 1:nSim){
    message('Starting iteration ', s)

    ## List of simulated panel datasets (one per cell)
    datDraw <- mclapply(paramList, simulate_one_cell, mc.cores=10)

    ## TODO Could run get_em_and_min_dist_estimates_random_initialization with use_md_as_initial_values_for_em=TRUE
    params0_hat <- mclapply(datDraw,
                            function(x) get_expectation_maximization_estimates(x, params0, max_iter=20, epsilon=0.001), mc.cores=10)

    estimated_deforestation_rates_ml <- sapply(params0_hat, function(p_hat) p_hat$P_list[[3]][1, 2])
    cor(estimated_deforestation_rates_ml, trueDeforestPr)
    mean(estimated_deforestation_rates_ml - trueDeforestPr)
                                                     
    viterbiList <- lapply(seq_along(datDraw),
                          function(i) apply_viterbi_path_in_parallel(datDraw[[i]],
                                                                     params_hat=params0_hat[[i]],
                                                                     max_cores=10))

    for(cell_idx in seq_along(datDraw)) {
        for(pixel_idx in seq_along(datDraw[[cell_idx]])) {
            datDraw[[cell_idx]][[pixel_idx]]$y_viterbi <- viterbiList[[cell_idx]][[pixel_idx]]
        }
    }

    ## We're using the EM/ML parameters to run viterbi (we could also try MD+viterbi)
    estimated_deforestation_rates_viterbi_ml <- sapply(datDraw, get_estimated_deforestation_rate_viterbi)
    estimated_deforestation_rates_observations <- sapply(datDraw, get_estimated_deforestation_rate_observations)
    estimated_deforestation_rates_ground_truth <- sapply(datDraw, get_estimated_deforestation_rate_ground_truth)

    ## Run regressions with various left hand side deforestation rates
    coef_observations <- coefficients(lm(estimated_deforestation_rates_observations ~ xVec))
    coef_ground_truth <- coefficients(lm(estimated_deforestation_rates_ground_truth ~ xVec))
    coef_viterbi_ml <- coefficients(lm(estimated_deforestation_rates_viterbi_ml ~ xVec))
    coef_ml <- coefficients(lm(estimated_deforestation_rates_ml ~ xVec))

    df <- data.frame(rbind(coef_observations, coef_viterbi_ml, coef_ml, coef_ground_truth))
    df$lhs_var <- str_replace(row.names(df), "coef_", "")

    names(df)[1] <- "alpha_hat"
    names(df)[2] <- "beta_hat"

    regResList[[s]] <- df
}

allRegDat <- rbindlist(regResList)
allRegDat <- melt(allRegDat, id.var="lhs_var")

filename <- sprintf("regression_simulation_estimated_coefficients_%s_sims.csv", nSim)
write.csv(allRegDat, filename, row.names=FALSE)

## Nicer order along x-axis
allRegDat$lhs_var <- factor(allRegDat$lhs_var, levels=c("observations", "viterbi_ml", "ml", "ground_truth"))

## Main MC graph
p <- ggplot(allRegDat, aes(x = lhs_var, y = value)) + geom_boxplot() + facet_wrap(~ variable) + theme_bw() + ylab("") + xlab("")

filename <- sprintf("regression_simulation_%s_sims.png", nSim)
ggsave(filename, plot=p, width=6, height=4, units="in")
