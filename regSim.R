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

cl <- makeCluster(8)

nCells <- 50  # TODO Bump back up to 100

nPixelPerCell <- 1000
alpha <- -3
beta <- 1

nSim <- 25  # TODO Bump back up to 100

set.seed(321321)

params0 <- get_params0()

xVec <- rnorm(nCells)  # Predictors that vary at the cell level
trueDeforestPr <- exp(alpha + beta * xVec) / (1+exp(alpha + beta * xVec))

tMax <- length(params0$P_list) + 1  # Number of time periods

## Deforestation only varies (with respect to xVec) in the 3rd period
## These parameters are held fixed across simulations
paramList <- list()
for (i in 1:length(xVec)) {
    paramList[[i]] <- params0
    paramList[[i]]$P_list[[3]][1,2] <- trueDeforestPr[i]
    paramList[[i]]$P_list[[3]][1,1] <- 1 - trueDeforestPr[i]
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

    stopifnot(all(c("point_id", "time", "y_viterbi") %in% names(dtable)))

    dtable[, y_viterbi_one_period_ahead := c(tail(y_viterbi, .N-1), NA), by="point_id"]

    ## Joint distribution of (Y_{t+1}, Y_{t})
    M_Y_joint_hat_list <- lapply(seq_len(max(dtable$time) - 1), function(fixed_t) {
        with(subset(dtable, time == fixed_t), prop.table(table(y_viterbi_one_period_ahead, y)))
    })

    return(get_transition_probs_from_M_S_joint(M_Y_joint_hat_list[[3]])[1, 2])
}

get_estimated_deforestation_rate_observations <- function(panel) {
    dtable <- rbindlist(Map(data.frame, panel))
    setkey(dtable, point_id)

    stopifnot(all(c("point_id", "time", "y") %in% names(dtable)))

    dtable[, y_one_period_ahead := c(tail(y, .N-1), NA), by="point_id"]

    ## Joint distribution of (Y_{t+1}, Y_{t})
    M_Y_joint_hat_list <- lapply(seq_len(max(dtable$time) - 1), function(fixed_t) {
        with(subset(dtable, time == fixed_t), prop.table(table(y_one_period_ahead, y)))
    })

    return(get_transition_probs_from_M_S_joint(M_Y_joint_hat_list[[3]])[1, 2])
}

regResList <- list()
for (s in 1:nSim){
    message('Starting iteration ', s)

    ## List of simulated panel datasets (one per cell)
    datDraw <- mclapply(paramList, simulate_one_cell, mc.cores=8)

    ## TODO Could run get_em_and_min_dist_estimates_random_initialization with use_md_as_initial_values_for_em=TRUE
    params0_hat <- mclapply(datDraw,
                            function(x) get_expectation_maximization_estimates(x, params0, max_iter=20, epsilon=0.001), mc.cores=8)

    estimated_deforestation_rates_ml <- sapply(params0_hat, function(p_hat) p_hat$P_list[[3]][1, 2])
    cor(estimated_deforestation_rates_ml, trueDeforestPr)
    mean(estimated_deforestation_rates_ml - trueDeforestPr)
                                                     
    viterbiList <- lapply(seq_along(datDraw),
                          function(i) apply_viterbi_path_in_parallel(datDraw[[i]],
                                                                     params_hat=params0_hat[[i]],
                                                                     max_cores=8))

    for(cell_idx in seq_along(datDraw)) {
        for(pixel_idx in seq_along(datDraw[[cell_idx]])) {
            datDraw[[cell_idx]][[pixel_idx]]$y_viterbi <- viterbiList[[cell_idx]][[pixel_idx]]
        }
    }

    ## Note: we could also try MD+viterbi
    estimated_deforestation_rates_viterbi_ml <- sapply(datDraw, get_estimated_deforestation_rate_viterbi)
    estimated_deforestation_rates_observations <- sapply(datDraw, get_estimated_deforestation_rate_observations)

    # TODO Also estimate deforestation rates using true hidden state

    mean(estimated_deforestation_rates_observations - trueDeforestPr)
    mean(estimated_deforestation_rates_viterbi_ml - trueDeforestPr)

    glm(estimated_deforestation_rates_observations ~ xVec, family = 'binomial')

    ## Run regressions with various left hand side deforestation rates
    coef_observations <- coefficients(glm(estimated_deforestation_rates_observations ~ xVec, family = 'binomial'))
    coef_viterbi_ml <- coefficients(glm(estimated_deforestation_rates_viterbi_ml ~ xVec, family = 'binomial'))
    coef_ml <- coefficients(glm(estimated_deforestation_rates_ml ~ xVec, family = 'binomial'))

    df <- data.frame(rbind(coef_observations, coef_viterbi_ml, coef_ml))
    df$lhs_var <- str_replace(row.names(df), "coef_", "")

    names(df)[1] <- "alpha_hat"
    names(df)[2] <- "beta_hat"

    regResList[[s]] <- df
}

allRegDat <- rbindlist(regResList)
allRegDat <- melt(allRegDat, id.var="lhs_var")

## Nicer order along x-axis
allRegDat$lhs_var <- factor(allRegDat$lhs_var, levels=c("observations", "viterbi_ml", "ml"))

## Main MC graph
title <- sprintf("Estimated Coefficients, %s Simulations\nEach simulation has %s cells, %s pixels per cell", nSim, nCells, nPixelPerCell)
p <- ggplot(allRegDat, aes(x = lhs_var, y = value)) + geom_boxplot() + facet_wrap(~ variable) + ggtitle(title)

filename <- sprintf("regression_simulation_%s_sims.png", nSim)
ggsave(filename, plot=p, width=6, height=4, units="in")
