rm(list = ls())
library(data.table)
library(optparse)
library(parallel)
source("hmm_functions.R")

## TODO There's a second R script that runs after this one, make sure it works, rename it, put it in readme

opt_list <- list(make_option("--n_simulations", default=100, type="integer"))
opt <- parse_args(OptionParser(option_list=opt_list))

message("command line options: ", paste(sprintf("%s=%s", names(opt), opt), collapse=", "))

simulation_df <- data.frame(n_points=c(100, 500, rep(1000, 13), 10000),
                            n_simulations=opt$n_simulations,
                            n_time_periods=c(rep(4, 10), 5, 6, rep(4,3), 4),
                            n_components = c(rep(2, 16)),
                            mu1 = c(rep(90, 16)), ## Put in as percent
                            defRt1 = c(rep(4, 16)), ## In as percent
                            defRtMid = c(rep(10, 16)), ## Percent
                            defRtLast = c(20, 20, seq(5, 40, by=5), rep(20, 6)),
                            prY11 = c(rep(90, 12), 75, 80, 95, 90), ## Correct classification Probability 1,1 in matrix
                            prY22 = c(rep(80, 16))) ## Correct classification Probability 2,2 in matrix

max_cores <- 4
set.seed(998877)

## Set output files for simulation
outfile_format <- paste0("baseline_simulation_", Sys.time(), "_iter_%s.rds")

## Output simulation parameters
iter_desc_outfile <- paste0("baseline_simulation_", Sys.time(), "_Desc.csv")
fwrite(simulation_df, file=iter_desc_outfile)

get_P_list <- function(deforestation_rates) {

    ## Returns a list of hidden state transition probabilities
    P_list <- lapply(deforestation_rates, function(deforestation_rate) {

        ## Note: the reforestation rate (second row of the transition probability matrix) is fixed
        P <- rbind(c(1 - deforestation_rate, deforestation_rate),
                   c(0.02, 0.98))

        return(P)
    })

    return(P_list)
}



run_single_simulation <- function(county_id, n_time_periods, n_points,n_components, mu1, prY11,prY22, deforestation_rates) {

    message(" simulating county_id ", county_id)

    params <- list(n_components=n_components,
                   mu=c(mu1, 1-mu1))

    ## Observation probabilities (rows are hidden states, columns are Y)
    ## Note: for now, every county has the same misclassification probabilities
    ## Note that this is the transpose of upsilon in the paper.
    pr_y <- rbind(c(prY11, 1 - prY11),
                  c(1 - prY22, prY22))

    P_list <- get_P_list(deforestation_rates)

    params$P_list <- P_list
    params$pr_y <- pr_y

    county_simulation <- replicate(n_points, simulate_hmm(params), simplify=FALSE)

    estimates <- get_em_and_min_dist_estimates_random_initialization(params=params, panel=county_simulation)

    return(list(simulation=county_simulation,
                estimates=estimates,
                params=params,
                deforestation_rates=deforestation_rates,
                id=county_id))
}


num_cores <- min(detectCores(), max_cores)
cluster <- makeCluster(num_cores)  # Call stopCluster when done

for (i in seq_len(nrow(simulation_df))) {

    n_simulations <- simulation_df$n_simulations[i]
    n_points <- simulation_df$n_points[i]
    n_time_periods <- simulation_df$n_time_periods[i]
    n_components <- simulation_df$n_components[i]
    mu1 <- simulation_df$mu1[i]/100
    prY11  <- simulation_df$prY11[i]/100
    prY22  <- simulation_df$prY22[i]/100
    deforestation_rates  <- c(simulation_df$defRt1[i]/100,
                              rep(simulation_df$defRtMid[i]/100, n_time_periods-3),
                              simulation_df$defRtLast[i]/100)

    message("Running simulation with n_simulations=", n_simulations, " n_points=", n_points, " n_time_periods=", n_time_periods)

    ## Note: only some of these objects change on every loop (e.g. n_points), others are constant
    clusterExport(cluster, c("baum_welch",
                             "eq_function_minimum_distance",
                             "get_P_list",
                             "get_deforestation_prob_from_P",
                             "get_expectation_maximization_estimates",
                             "get_em_and_min_dist_estimates_random_initialization",
                             "get_min_distance_estimates",
                             "get_random_initial_parameters",
                             "get_transition_probs_from_M_S_joint",
                             "is_diag_dominant",
                             "mu1",
                             "n_components",
                             "n_points",
                             "n_time_periods",
                             "prY11",
                             "prY22",
                             "deforestation_rates",
                             "objfn_minimum_distance",
                             "simulate_discrete_markov",
                             "simulate_hmm",
                             "run_single_simulation",
                             "valid_panel_element",
                             "valid_parameters"))
    
    simulations <- parLapply(cluster, seq_len(n_simulations), function(n) {
        run_single_simulation(county_id=n, n_time_periods=n_time_periods, n_points=n_points,
                              n_components = n_components,mu1 = mu1, prY11 = prY11, prY22 = prY22,
                              deforestation_rates = deforestation_rates)
    })

    outfile <- sprintf(outfile_format, i)
    message("Saving ", outfile)
    saveRDS(simulations, file=outfile)
}

stopCluster(cluster)
