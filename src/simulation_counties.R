library(data.table)

source("hmm_functions.R")

n_time_periods <- 4
n_points_per_county <- 1000
n_counties <- 50

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

simulate_single_county <- function(county_id, n_time_periods, n_points_per_county) {

    county_fixed_effect <- rnorm(1)

    ## Note: X is an observed predictor that varies by (county, time)
    county_X <- rnorm(n_time_periods)

    ## Note: for now, every county has the same initial distribution over hidden states
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

    return(list(simulation=county_simulation,
                params=county_params,
                fixed_effect=county_fixed_effect,
                X=county_X,
                id=county_id))
}

counties <- lapply(seq_len(n_counties), function(n) {
    simulate_single_county(county_id=n, n_time_periods=n_time_periods, n_points_per_county=n_points_per_county)
})

county_dfs <- lapply(counties, function(county) {

    ## TODO Also return naive, MD and EM deforestation probability estimates

    true_deforestation_probability <- sapply(county$params$P_list, function(P) {
        return(P[1, 2])
    })
    return(data.table(x=county$X,
                      fixed_effect=county$fixed_effect,
                      time=seq_along(county$X),
                      county_id=county$id,
                      true_deforestation_probability=true_deforestation_probability))
})

county_df <- rbindlist(county_dfs)
setkey(county_df, county_id)

county_df[, county_id_factor := factor(county_id)]
county_df[, true_y := log(true_deforestation_probability)]

summary(lm(true_y ~ x + county_id_factor, data=county_df))

# Note: this should recover the coefficients in get_deforestation_probability
summary(lm(true_y ~ x + fixed_effect, data=county_df))

## TODO Estimate relationship between county_X and deforestation rate (using naive, EM, and MD estimators)
