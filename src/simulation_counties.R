rm(list = ls())
library(data.table)
library(stargazer)
library(parallel)
source("hmm_functions.R")


## '2019-09-24 19:24:42'
## simulation_df <- data.frame(n_points_per_county=c(100, 500, 1000, rep(1000,8)),
##                             n_counties=c(50, 50, 50, rep(50,8)),
##                             n_time_periods=c(4, 4, 4 , rep(4, 8)),
##                             n_components = c(2, 2, 2, rep(2, 8)),
##                             mu1 = c(90, 90, 90, rep(90,8)), ##Put in as percent
##                             defRt1 = c(4, 4, 4, rep(4,8)), ##In as percent
##                             defRtMid = c(10, 10, 10, rep(10,8)), ##same
##                             defRtLast = c(20, 20, 20, seq(5,40,by=5)),
##                             prY11 = c(90, 90, 90,rep(90,8)), ##Correct classification Probability 1,1 in matrix
##                             prY22 = c(80, 80, 80, rep(80,8))) ##Correct classification Probability 2,2 in matrix

simulation_df <- data.frame(n_points_per_county=c( rep(1000,6)),
                            n_counties=c(rep(50,6)),
                            n_time_periods=c(5, 6, rep(4,4)),
                            n_components = c(rep(2,6)),
                            mu1 = c(rep(90,6)), ##Put in as percent
                            defRt1 = c(rep(4,6)), ##In as percent
                            defRtMid = c(rep(10,6)), ##same
                            defRtLast = c(rep(20,6)),
                            prY11 = c(rep(90,2),75,80,85,95), ##Misclassification Probability 1,1 in matrix
                            prY22 = c(rep(80,6))) ##Misclassification Probability 2,2 in matrix


max_cores <- 4
set.seed(998877)


##Generate matrix of simulation values (constant across all simulations)
maxNCount  <- max(simulation_df$n_counties)
maxT <- max(simulation_df$n_time_periods)
feVec  <- vector(mode='numeric',length = maxNCount)*0

##Set output files for simulation
county_outfile_format <- paste0("county_simulation_",Sys.time(),"_iter_%s.rds")

##Output simulation parameters
iter_desc_outfile <- paste0("county_simulation_",Sys.time(),"_Desc.csv")
fwrite(simulation_df,file = iter_desc_outfile)


## Keep track of the estimates from the county regression
simulation_df$county_regression_naive_beta_hat <- NA
simulation_df$county_regression_naive_beta_hat_std_err <- NA
simulation_df$county_regression_em_beta_hat <- NA
simulation_df$county_regression_em_beta_hat_std_err <- NA
simulation_df$county_regression_md_beta_hat <- NA
simulation_df$county_regression_md_beta_hat_std_err <- NA

#regression_summary_outfile_format <- "county_simulation_regressions_%s_points_%s_counties.tex"
#simulation_df_outfile <- "county_simulations_summary.tex"




get_deforestation_probability <- function(x, county_fixed_effect) {

    return(x)

}

# TODO Make sure the distribution of county_X and county_fixed_effect leads to somewhat realistic deforestation probabilities
#curve(get_deforestation_probability(x, county_fixed_effect=1), from=-2, to=2)
#curve(get_deforestation_probability(x, county_fixed_effect=0), from=-2, to=2, lty=2, col="red", add=TRUE)

get_P_list <- function(county_fixed_effect, county_X) {

    ## Returns a list of hidden state transition probabilities
    P_list <- lapply(county_X, function(x) {

        deforestation <- get_deforestation_probability(x, county_fixed_effect)

        ## Note: the reforestation rate is constant across time and counties for now
        ## (second row of the transition probability matrix is fixed)
        P <- rbind(c(1 - deforestation, deforestation),
                   c(0.02, 0.98))

        return(P)
    })

    return(P_list)
}



simulate_single_county <- function(county_id, n_time_periods, n_points_per_county,n_components, mu1, prY11,prY22, county_fixed_effect, county_X) {

    message(" simulating county_id ", county_id)


    ## Note: for now, every county has the same initial distribution over hidden states (mu)
    county_params <- list(n_components=n_components,
                          mu=c(mu1, 1-mu1))

    ## Observation probabilities (rows are hidden states, columns are Y)
    ## Note: for now, every county has the same misclassification probabilities
    pr_y <- rbind(c(prY11, 1-prY11),
                  c(1-prY22, prY22))

    P_list <- get_P_list(county_fixed_effect, county_X)

    county_params$P_list <- P_list
    county_params$pr_y <- pr_y

    county_simulation <- replicate(n_points_per_county, simulate_hmm(county_params), simplify=FALSE)

    estimates <- get_hmm_and_minimum_distance_estimates_random_initialization(params=county_params, panel=county_simulation)

    return(list(simulation=county_simulation,
                estimates=estimates,
                params=county_params,
                fixed_effect=county_fixed_effect,
                X=county_X,
                id=county_id))
}


num_cores <- min(detectCores(), max_cores)
cluster <- makeCluster(num_cores)  # Call stopCluster when done

for (i in seq_len(nrow(simulation_df))) {

    n_counties <- simulation_df$n_counties[i]
    n_points_per_county <- simulation_df$n_points_per_county[i]
    n_time_periods <- simulation_df$n_time_periods[i]
    n_components <- simulation_df$n_components[i]
    mu1 <- simulation_df$mu1[i]/100
    prY11  <- simulation_df$prY11[i]/100
    prY22  <- simulation_df$prY22[i]/100
    xVec  <- c(simulation_df$defRt1[i]/100,
               rep(simulation_df$defRtMid[i]/100,n_time_periods-3),
               simulation_df$defRtLast[i]/100)
               

    message("Running simulation with n_counties=", n_counties, " n_points_per_county=", n_points_per_county, " n_time_periods=", n_time_periods)

    ## Note: only some of these objects change on every loop (e.g. n_points_per_county), others are constant
    clusterExport(cluster, c("baum_welch",
                             "eq_function_minimum_distance",
                             "get_P_list",
                             "get_deforestation_prob_from_P",
                             "get_deforestation_probability",
                             "get_expectation_minimization_estimates",
                             "get_hmm_and_minimum_distance_estimates_random_initialization",
                             "get_min_distance_estimates",
                             "get_random_initial_parameters",
                             "get_transition_probs_from_M_S_joint",
                             "n_points_per_county",
                             "n_time_periods",
                             "n_components",
                             "mu1",
                             "prY11",
                             "prY22",
                             "feVec",
                             "xVec",
                             "objfn_minimum_distance",
                             "simulate_discrete_markov",
                             "simulate_hmm",
                             "simulate_single_county",
                             "valid_panel_element",
                             "valid_parameters"))


    counties <- parLapply(cluster, seq_len(n_counties), function(n) {
        simulate_single_county(county_id=n, n_time_periods=n_time_periods, n_points_per_county=n_points_per_county,
                               n_components = n_components,mu1 = mu1, prY11 = prY11, prY22 = prY22,
                               county_fixed_effect = feVec[n], county_X = xVec)
    })

    ## ##Debug
    ## minDistMiscPr <- sapply(counties, function(x) x$estimates$min_dist_params_hat_best_objfn$pr_y[1,1])    
    ## emMiscPr  <- sapply(counties, function(x) x$estimates$em_params_hat_best_likelihood$pr_y[1,1])
    ## trueMiscPr  <- sapply(counties, function(x) x$params$pr_y[1,1])

    ## mdErrMisc  <- mean(abs(minDistMiscPr-trueMiscPr))
    ## emErrMisc  <- mean(abs(emMiscPr-trueMiscPr))
    
    ## mdErrMat  <- matrix(nrow=5,ncol=2)
    ## emErrMat <- matrix(nrow = 5, ncol = 2)
    ## mdWorseEmMat <- matrix(nrow = 5, ncol = 2)
    ## for (ti in 1:5){
    ##     for (ro in 1:2){
    ##         minDistTransPr <- sapply(counties, function(x) x$estimates$min_dist_params_hat_best_objfn$P_list[[ti]][ro,1])    
    ##         emTransPr  <- sapply(counties, function(x) x$estimates$em_params_hat_best_likelihood$P_list[[ti]][ro,1])
    ##         trueTransPr  <- sapply(counties, function(x) x$params$P_list[[ti]][ro,1])
    ##         mdErrMat[ti,ro]  <- mean(abs(minDistTransPr - trueTransPr))
    ##         emErrMat[ti,ro]  <- mean(abs(emTransPr - trueTransPr))
    ##         mdWorseEmMat[ti,ro]  <- mean(abs(minDistTransPr - trueTransPr) -
    ##                                      abs(emTransPr - trueTransPr) > .01)
    ##     }
    ## }
    

    
    ## county_dfs <- lapply(counties, get_data_table_summarizing_single_county_simulation)

    ## county_df <- rbindlist(county_dfs)

    ## ## Note: this sorts county_df first by county_id and then by time,
    ## ## the sort order is necessary for computing first differences correctly
    ## setkey(county_df, county_id, time)

    ## county_df[, county_id_factor := factor(county_id)]
    ## county_df[, y_true := log(true_deforestation_probability / (1 - true_deforestation_probability))]
    ## county_df[, y_naive:= log(estimated_deforestation_probability_naive / (1 - estimated_deforestation_probability_naive))]
    ## county_df[, y_em := log(estimated_deforestation_probability_em / (1 - estimated_deforestation_probability_em))]
    ## county_df[, y_md := log(estimated_deforestation_probability_md / (1 - estimated_deforestation_probability_md))]

    ## county_df[, x_first_diff := c(NA, diff(x)), by="county_id"]
    ## county_df[, y_true_first_diff := c(NA, diff(y_true)), by="county_id"]
    ## county_df[, y_naive_first_diff := c(NA, diff(y_naive)), by="county_id"]
    ## county_df[, y_em_first_diff := c(NA, diff(y_em)), by="county_id"]
    ## county_df[, y_md_first_diff := c(NA, diff(y_md)), by="county_id"]

    ## Note: if you want to skip the simulation,
    ## you can load county_df_outfile and run the regressions (lm) in the lines below
    county_outfile <- sprintf(county_outfile_format, i)
    message("Saving ", county_outfile)
    saveRDS(counties, file = county_outfile)
}


stopCluster(cluster)
