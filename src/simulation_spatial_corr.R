library(data.table)
library(ggplot2)

source("hmm_functions.R")
source("hmm_parameters.R")

set.seed(321123)

params0 <- get_params0()
params1 <- get_params1(params0)

## In this simulation, we want the hidden state S_{it}
## to vary at the _field_ level rather than the pixel level
## However, classification errors will occur at the pixel level,
## and we will use pixel-level observations to estimate the model's parameters
n_fields <- 100

fields <- lapply(seq_len(n_fields), function(field_id, params=params0) {
    state <- simulate_discrete_markov(params)
    return(list(state=state, field_id=field_id))
})

## At 30 meter resolution, this would mean 30^2 * 200^2 meters^2 of area,
## or roughly 36 kilometers^2 (8,896 acres)
n_pixels_per_side <- 200

## Construct a join table for going from pixels to their field_id
df <- expand.grid(pixel_i=seq_len(n_pixels_per_side), pixel_j=seq_len(n_pixels_per_side),
                  field_id=0)

n_fields_per_side <- sqrt(n_fields)

## TODO Could have fields of different sizes by changing the edges here
field_edges <- expand.grid(pixel_i=seq(1, n_pixels_per_side, n_pixels_per_side / n_fields_per_side),
                           pixel_j=seq(1, n_pixels_per_side, n_pixels_per_side / n_fields_per_side))
stopifnot(nrow(field_edges) == n_fields)

for(field_id in seq_len(n_fields)) {
    cutoff_i <- field_edges[field_id, ]$pixel_i
    cutoff_j <- field_edges[field_id, ]$pixel_j
    df[df$pixel_i >= cutoff_i & df$pixel_j >= cutoff_j, ]$field_id <- field_id
}

## How many pixels are in each field?  Should be (n_pixels_per_side / n_fields_per_side) ^ 2 pixels per field
## TODO Could generalize and have fields of different sizes
table(df$field_id)

pixel_panel <- lapply(seq_len(nrow(df)), function(row_index, params=params0) {
    field_id <- df[row_index, ]$field_id
    pixel_i <- df[row_index, ]$pixel_i
    pixel_j <- df[row_index, ]$pixel_j
    field_state <- fields[[field_id]]$state
    y <- vapply(field_state, function(s) {
        sample(seq_len(ncol(params$pr_y)), size=1, prob=params$pr_y[s, ])
    }, FUN.VALUE=1)
    time <- seq_along(y)
    return(list(field_state=field_state,
                field_id=field_id,
                y=y,
                time=time,
                pixel_i=pixel_i,
                pixel_j=pixel_j))
})

## Inspect the observations y and hidden state at the first two pixels
pixel_panel[[1]]
pixel_panel[[2]]

## Initialize EM at true parameter values (easy case)
params0_hat <- get_expectation_maximization_estimates(pixel_panel, params0, max_iter=20, epsilon=0.001)

stopifnot(all(diff(params0_hat$loglik) > 0))  # Loglik should be increasing
max(abs(c(params0_hat$P_list, recursive=TRUE) - c(params0$P_list, recursive=TRUE)))  # Largest error in time-varying transition probabilities
max(abs(params0_hat$pr_y - params0$pr_y))  # Largest error in observation probabilities (aka misclassification probabilities)

## Initialize EM at incorrect parameter values (more difficult)
params1_hat <- get_expectation_maximization_estimates(pixel_panel, params1, max_iter=20, epsilon=0.001)

max(abs(c(params1_hat$P_list, recursive=TRUE) - c(params0$P_list, recursive=TRUE)))  # Largest error in time-varying transition probabilities
max(abs(params1_hat$pr_y - params0$pr_y))  # Largest error in observation probabilities (aka misclassification probabilities)

dtable <- rbindlist(Map(data.frame, pixel_panel))
head(dtable)

dtable[, time_label := sprintf("time %s", time)]

colors <- c("#dfc27d", "#018571")

p <- (ggplot(subset(dtable, time == 1), aes(x=pixel_i, y=pixel_j, fill=factor(field_state))) +
      geom_raster() +
      scale_fill_manual("true field state", values=colors) +
      xlab("pixel coordinate (easting)") +
      ylab("pixel coordinate (northing)"))
p

p <- (ggplot(dtable, aes(x=pixel_i, y=pixel_j, fill=factor(field_state))) +
      facet_wrap(~ time_label) +
      geom_raster() +
      scale_fill_manual("true field state", values=colors) +
      xlab("pixel coordinate (easting)") +
      ylab("pixel coordinate (northing)"))
p
ggsave("simulation_spatial_corr_true_field_state.png", plot=p, width=6, height=4, units="in")

p <- (ggplot(subset(dtable, time == 1), aes(x=pixel_i, y=pixel_j, fill=factor(y))) +
      geom_raster() +
      scale_fill_manual("observed y", values=colors) +
      xlab("pixel coordinate (easting)") +
      ylab("pixel coordinate (northing)"))
p

p <- (ggplot(dtable, aes(x=pixel_i, y=pixel_j, fill=factor(y))) +
      facet_wrap(~ time_label) +
      geom_raster() +
      scale_fill_manual("observed y", values=colors) +
      xlab("pixel coordinate (easting)") +
      ylab("pixel coordinate (northing)"))
p
ggsave("simulation_spatial_corr_observed_y.png", plot=p, width=6, height=4, units="in")

dtable[, pixel_id := sprintf("pixel_%s_%s", pixel_i, pixel_j)]

dtable[, y_one_period_ahead := c(tail(y, .N-1), NA), by="pixel_id"]
dtable[, y_two_periods_ahead := c(tail(y, .N-2), NA, NA), by="pixel_id"]
head(dtable[, c("pixel_id", "time", "y", "y_one_period_ahead", "y_two_periods_ahead"), with=FALSE], 25)  # Sanity check

## Compute the joint distribution of (Y_{t+1}, Y_{t}) for each time period t
M_Y_joint_hat <- lapply(seq_len(max(dtable$time) - 1), function(fixed_t) {
    with(subset(dtable, time == fixed_t), prop.table(table(y_one_period_ahead, y)))
})

## Naive estimate of transitions using observed Y
P1_hat_naive <- get_transition_probs_from_M_S_joint(M_Y_joint_hat[[1]])
P2_hat_naive <- get_transition_probs_from_M_S_joint(M_Y_joint_hat[[2]])
P3_hat_naive <- get_transition_probs_from_M_S_joint(M_Y_joint_hat[[3]])

## Naive estimates have much larger errors than EM estimates
max(abs(P1_hat_naive - params0$P_list[[1]]))
max(abs(params0_hat$P_list[[1]] - params0$P_list[[1]]))

max(abs(P2_hat_naive - params0$P_list[[2]]))
max(abs(params0_hat$P_list[[2]] - params0$P_list[[2]]))

max(abs(P3_hat_naive - params0$P_list[[3]]))
max(abs(params0_hat$P_list[[3]] - params0$P_list[[3]]))
