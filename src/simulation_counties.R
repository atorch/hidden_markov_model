library(data.table)
library(stargazer)
library(parallel)

source("hmm_functions.R")


## Modify this dataframe to vary the number of simulations and the simulation parameters
simulation_df <- data.frame(n_points_per_county=c(1000, 2000, 4000),
                            n_counties=c(200, 100, 50),
                            n_time_periods=c(10, 10, 10))

## Keep track of the estimates from the county regression
simulation_df$county_regression_naive_beta_hat <- NA
simulation_df$county_regression_naive_beta_hat_std_err <- NA
simulation_df$county_regression_em_beta_hat <- NA
simulation_df$county_regression_em_beta_hat_std_err <- NA
simulation_df$county_regression_md_beta_hat <- NA
simulation_df$county_regression_md_beta_hat_std_err <- NA

county_df_outfile_format <- "county_simulation_%s_points_%s_counties.csv"
regression_summary_outfile_format <- "county_simulation_regressions_%s_points_%s_counties.tex"
simulation_df_outfile <- "county_simulations_summary.tex"

max_cores <- 12

set.seed(998877)

get_deforestation_probability <- function(x, county_fixed_effect) {

    ## Note: x is a county- and time-varying predictor
    return(1 / (1 + exp(-(-4 + x + county_fixed_effect))))

}

# TODO Make sure the distribution of county_X and county_fixed_effect leads to somewhat realistic deforestation probabilities
curve(get_deforestation_probability(x, county_fixed_effect=1), from=-2, to=2)
curve(get_deforestation_probability(x, county_fixed_effect=0), from=-2, to=2, lty=2, col="red", add=TRUE)

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

get_deforestation_prob_from_P <- function(P) {
    ## Note: this assumes we have 2 hidden states, and state 1 is forest
    return(P[1, 2])
}

simulate_single_county <- function(county_id, n_time_periods, n_points_per_county) {

    message(" simulating county_id ", county_id)

    county_fixed_effect <- rnorm(1)

    ## Note: X is an observed predictor that varies by (county, time)
    county_X <- rnorm(n_time_periods)

    ## Note: for now, every county has the same initial distribution over hidden states (mu)
    county_params <- list(n_components=2,
                          mu=c(0.9, 0.1))

    ## Observation probabilities (rows are hidden states, columns are Y)
    ## Note: for now, every county has the same misclassification probabilities
    pr_y <- rbind(c(0.90, 0.10),
                  c(0.20, 0.80))

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

get_data_table_summarizing_single_county_simulation <- function(county) {

    true_deforestation_prob <- sapply(county$params$P_list,
                                      get_deforestation_prob_from_P)

    estimated_deforestation_prob_em <- sapply(county$estimates$em_params_hat_best_likelihood$P_list,
                                              get_deforestation_prob_from_P)

    estimated_deforestation_prob_md <- sapply(county$estimates$min_dist_params_hat_best_objfn$P_list,
                                              get_deforestation_prob_from_P)

    P_hat_naive <- lapply(county$estimates$M_Y_joint_hat, get_transition_probs_from_M_S_joint)
    estimate_deforestation_prob_naive <- sapply(P_hat_naive, get_deforestation_prob_from_P)

    ## EM means expectation-maximization, MD means minimum distance
    return(data.table(x=county$X,
                      fixed_effect=county$fixed_effect,
                      time=seq_along(county$X),
                      county_id=county$id,
                      true_deforestation_probability=true_deforestation_prob,
                      estimated_deforestation_probability_naive=estimate_deforestation_prob_naive,
                      estimated_deforestation_probability_em=estimated_deforestation_prob_em,
                      estimated_deforestation_probability_md=estimated_deforestation_prob_md))
}

num_cores <- min(detectCores(), max_cores)
cluster <- makeCluster(num_cores)  # Call stopCluster when done

for (i in seq_len(nrow(simulation_df))) {

    n_counties <- simulation_df$n_counties[i]
    n_points_per_county <- simulation_df$n_points_per_county[i]
    n_time_periods <- simulation_df$n_time_periods[i]

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
                             "objfn_minimum_distance",
                             "simulate_discrete_markov",
                             "simulate_hmm",
                             "simulate_single_county",
                             "valid_panel_element",
                             "valid_parameters"))

    counties <- parLapply(cluster, seq_len(n_counties), function(n) {
        simulate_single_county(county_id=n, n_time_periods=n_time_periods, n_points_per_county=n_points_per_county)
    })

    county_dfs <- lapply(counties, get_data_table_summarizing_single_county_simulation)

    county_df <- rbindlist(county_dfs)

    ## Note: this sorts county_df first by county_id and then by time,
    ## the sort order is necessary for computing first differences correctly
    setkey(county_df, county_id, time)

    county_df[, county_id_factor := factor(county_id)]
    county_df[, y_true := log(true_deforestation_probability / (1 - true_deforestation_probability))]
    county_df[, y_naive:= log(estimated_deforestation_probability_naive / (1 - estimated_deforestation_probability_naive))]
    county_df[, y_em := log(estimated_deforestation_probability_em / (1 - estimated_deforestation_probability_em))]
    county_df[, y_md := log(estimated_deforestation_probability_md / (1 - estimated_deforestation_probability_md))]

    county_df[, x_first_diff := c(NA, diff(x)), by="county_id"]
    county_df[, y_true_first_diff := c(NA, diff(y_true)), by="county_id"]
    county_df[, y_naive_first_diff := c(NA, diff(y_naive)), by="county_id"]
    county_df[, y_em_first_diff := c(NA, diff(y_em)), by="county_id"]
    county_df[, y_md_first_diff := c(NA, diff(y_md)), by="county_id"]

    ## Note: if you want to skip the simulation,
    ## you can load county_df_outfile and run the regressions (lm) in the lines below
    county_df_outfile <- sprintf(county_df_outfile_format, n_points_per_county, n_counties)
    message("Saving ", county_df_outfile)
    write.csv(county_df, county_df_outfile, row.names=FALSE)

    with(county_df, plot(true_deforestation_probability, estimated_deforestation_probability_naive)); abline(a=0, b=1, lty=2)
    with(county_df, plot(true_deforestation_probability, estimated_deforestation_probability_em)); abline(a=0, b=1, lty=2)
    with(county_df, plot(true_deforestation_probability, estimated_deforestation_probability_md)); abline(a=0, b=1, lty=2)

    ## Note: the constant term (intercept) in the county_id_factor regression includes the estimated fixed effect for county 1
    ## This should (exactly) recover the coefficients in get_deforestation_probability, hence the
    ## "essentially perfect fit: summary may be unreliable" warning
    ## Note: the 0 in the first diff regression formula means there is no constant (since there are no time effects in get_deforestation_probability)
    summary(lm(y_true ~ x + county_id_factor, data=county_df))
    summary(lm(y_true ~ x + fixed_effect, data=county_df))
    summary(lm(y_true_first_diff ~ 0 + x_first_diff, data=county_df))

    ## Note: the coefficients on both x and fixed_effect appear biased relative to the true coefficients in get_deforestation_probability
    summary(lm(y_naive ~ x + fixed_effect, data=county_df))
    summary(lm(y_naive ~ x + county_id_factor, data=county_df))
    summary(lm(y_naive_first_diff ~ 0 + x_first_diff, data=county_df))

    ## Looks like both y_em and y_md work better than y_naive,
    ## but why are the std errors so much larger (especially for y_md)?
    summary(lm(y_em ~ x + fixed_effect, data=county_df))
    summary(lm(y_em ~ x + county_id_factor, data=county_df))
    summary(lm(y_em_first_diff ~ 0 + x_first_diff, data=county_df))

    summary(lm(y_md ~ x + fixed_effect, data=county_df))
    summary(lm(y_md ~ x + county_id_factor, data=county_df))
    summary(lm(y_md_first_diff ~ 0 + x_first_diff, data=county_df))

    ## Save model objects so that we can pass them to stargazer (and save latex tables)
    model_naive <- lm(y_naive ~ x + county_id_factor, data=county_df)
    model_em <- lm(y_em ~ x + county_id_factor, data=county_df)
    model_md <- lm(y_md ~ x + county_id_factor, data=county_df)

    regression_summary_outfile <- sprintf(regression_summary_outfile_format, n_points_per_county, n_counties)
    title <- sprintf("Regressions with %s points per county, %s counties", n_points_per_county, n_counties)
    message("Saving ", regression_summary_outfile)
    stargazer(model_naive,
              model_em,
              model_md,
              keep=c("x"),
              title=title,
              out=regression_summary_outfile,
              notes=c("Estimates of county fixed effects are omitted"))

    coef_naive <- coef(summary(model_naive))
    coef_em <- coef(summary(model_em))
    coef_md <- coef(summary(model_md))

    simulation_df$county_regression_naive_beta_hat[i] <- coef_naive["x", "Estimate"]
    simulation_df$county_regression_naive_beta_hat_std_err[i] <- coef_naive["x", "Std. Error"]

    simulation_df$county_regression_em_beta_hat[i] <- coef_em["x", "Estimate"]
    simulation_df$county_regression_em_beta_hat_std_err[i] <- coef_em["x", "Std. Error"]

    simulation_df$county_regression_md_beta_hat[i] <- coef_md["x", "Estimate"]
    simulation_df$county_regression_md_beta_hat_std_err[i] <- coef_md["x", "Std. Error"]

}

stopCluster(cluster)

stargazer(simulation_df,
          summary=FALSE,
          rownames=FALSE,
          out=simulation_df_outfile)
