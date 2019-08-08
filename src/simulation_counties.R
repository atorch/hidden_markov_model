library(data.table)

source("hmm_functions.R")

n_time_periods <- 10
n_points_per_county <- 2000
n_counties <- 20

county_df_outfile <- "county_simulation.csv"

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

counties <- lapply(seq_len(n_counties), function(n) {
    simulate_single_county(county_id=n, n_time_periods=n_time_periods, n_points_per_county=n_points_per_county)
})

county_dfs <- lapply(counties, function(county) {

    true_deforestation_prob <- sapply(county$params$P_list,
                                      get_deforestation_prob_from_P)

    estimated_deforestation_prob_hmm <- sapply(county$estimates$hmm_params_hat_best_likelihood$P_list,
                                               get_deforestation_prob_from_P)

    estimated_deforestation_prob_md <- sapply(county$estimates$min_dist_params_hat_best_objfn$P_list,
                                              get_deforestation_prob_from_P)

    P_hat_naive <- lapply(county$estimates$M_Y_joint_hat, get_transition_probs_from_M_S_joint)
    estimate_deforestation_prob_naive <- sapply(P_hat_naive, get_deforestation_prob_from_P)

    return(data.table(x=county$X,
                      fixed_effect=county$fixed_effect,
                      time=seq_along(county$X),
                      county_id=county$id,
                      true_deforestation_probability=true_deforestation_prob,
                      estimated_deforestation_probability_naive=estimate_deforestation_prob_naive,
                      estimated_deforestation_probability_hmm=estimated_deforestation_prob_hmm,
                      estimated_deforestation_probability_md=estimated_deforestation_prob_md))
})

county_df <- rbindlist(county_dfs)
setkey(county_df, county_id)

county_df[, county_id_factor := factor(county_id)]
county_df[, y_true := log(true_deforestation_probability / (1 - true_deforestation_probability))]
county_df[, y_naive:= log(estimated_deforestation_probability_naive / (1 - estimated_deforestation_probability_naive))]
county_df[, y_hmm := log(estimated_deforestation_probability_hmm / (1 - estimated_deforestation_probability_hmm))]
county_df[, y_md := log(estimated_deforestation_probability_md / (1 - estimated_deforestation_probability_md))]

## Note: if you want to skip the simulation,
## you can load county_df_outfile and run the regressions (lm) in the lines below
write.csv(county_df, county_df_outfile, row.names=FALSE)

with(county_df, plot(true_deforestation_probability, estimated_deforestation_probability_naive)); abline(a=0, b=1, lty=2)
with(county_df, plot(true_deforestation_probability, estimated_deforestation_probability_hmm)); abline(a=0, b=1, lty=2)
with(county_df, plot(true_deforestation_probability, estimated_deforestation_probability_md)); abline(a=0, b=1, lty=2)

## Note: the constant term (intercept) in this regression includes the estimated fixed effect for county 1
summary(lm(y_true ~ x + county_id_factor, data=county_df))

## Note: this should (exactly) recover the coefficients in get_deforestation_probability
summary(lm(y_true ~ x + fixed_effect, data=county_df))

## Note: the coefficients on both x and fixed_effect appear biased relative to the true coefficients in get_deforestation_probability
summary(lm(y_naive ~ x + fixed_effect, data=county_df))

## TODO Do these work better than y_naive?
## Looks like they do, but why are the std errors so much larger (especially for y_md)?
summary(lm(y_hmm ~ x + fixed_effect, data=county_df))
summary(lm(y_md ~ x + fixed_effect, data=county_df))
