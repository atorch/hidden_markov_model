source("hmm_functions.R")

n_time_periods <- 4
n_points_per_county <- 1000
n_counties <- 200

get_deforestation_probability <- function(x, county_fixed_effect) {

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

    return(list(county_simulation=county_simulation,
                county_params=county_params,
                county_fixed_effect=county_fixed_effect,
                county_X=county_X))
}

simulate_single_county(county_id=0, n_time_periods=n_time_periods, n_points_per_county=n_points_per_county)

## TODO Simulate multiple counties
## Estimate relationship between county_X and deforestation rate (using naive, EM, and MD estimators)
