library(data.table)
library(ggplot2)
library(parallel)

source("hmm_functions.R")
source("hmm_parameters.R")
source("ising.R")

set.seed(321123)

run_single_simulation <- function(simulation_id, params0, adjacency, n_pixels_per_side) {

    ## In this simulation, we want the hidden state S_{it}
    ## to vary at the _field_ level rather than the pixel level
    ## However, classification errors will occur at the pixel level,
    ## and we will use pixel-level observations to estimate the model's parameters

    n_fields <- params0$n_fields

    fields <- lapply(seq_len(n_fields), function(field_id, params=params0) {
        state <- simulate_discrete_markov(params)
        return(list(state=state, field_id=field_id))
    })

    ## Construct a join table for going from pixels to their field_id
    df <- expand.grid(pixel_i=seq_len(n_pixels_per_side),
                      pixel_j=seq_len(n_pixels_per_side),
                      field_id=0)

    n_fields_per_side <- sqrt(n_fields)

    field_edges <- expand.grid(pixel_i=seq(1, n_pixels_per_side, n_pixels_per_side / n_fields_per_side),
                               pixel_j=seq(1, n_pixels_per_side, n_pixels_per_side / n_fields_per_side))
    stopifnot(nrow(field_edges) == n_fields)

    for(field_id in seq_len(n_fields)) {
        cutoff_i <- field_edges[field_id, ]$pixel_i
        cutoff_j <- field_edges[field_id, ]$pixel_j
        df[df$pixel_i >= cutoff_i & df$pixel_j >= cutoff_j, ]$field_id <- field_id
    }

    if(params0$include_z_in_simulation) {
        if(params0$z_constant_over_time) {
            df$z <- simulate_ising(n_pixels=nrow(df), adjacency=adjacency, beta=params0$ising_beta, n_iter=50)
        } else {
            ## If Z is not constant over time, it is i.i.d. over time
            n_time_periods <- length(fields[[1]]$state)
            for(time in seq_len(n_time_periods)) {
                z_colname <- sprintf("z_%s", time)
                df[, z_colname] <- simulate_ising(n_pixels=nrow(df), adjacency=adjacency, beta=params0$ising_beta, n_iter=50)
            }
        }
    }

    pixel_panel <- lapply(seq_len(nrow(df)), function(row_index, params=params0) {

        field_id <- df[row_index, ]$field_id
        pixel_i <- df[row_index, ]$pixel_i
        pixel_j <- df[row_index, ]$pixel_j

        ## If Z is included in the simulation, construct a vector z
        ## of the same length as the hidden state (i.e. one value per time period)
        n_time_periods <- length(fields[[field_id]]$state)
        if(params$include_z_in_simulation) {
            if(params0$z_constant_over_time) {
                z <- as.numeric(rep(df[row_index, ]$z, n_time_periods))
            } else {
                z_colnames <- sapply(seq_len(n_time_periods), function(t) sprintf("z_%s", t))
                z <- as.numeric(df[row_index, z_colnames])
            }
        } else {
            z <- rep(NA, n_time_periods)
        }
        stopifnot(length(z) == n_time_periods)

        y <- vapply(seq_len(n_time_periods), function(time) {
            s <- fields[[field_id]]$state[time]
            if(params$include_z_in_simulation) {
                if(z[time] == 1) {
                    sample(seq_len(ncol(params$pr_y)), size=1, prob=params$pr_y_given_cloudy[s, ])
                } else {
                    sample(seq_len(ncol(params$pr_y)), size=1, prob=params$pr_y_given_clear[s, ])
                }
            } else  {
                sample(seq_len(ncol(params$pr_y)), size=1, prob=params$pr_y[s, ])
            }
        }, FUN.VALUE=1)
        time <- seq_along(y)

        if(params$pr_missing_data > 0) {
            ## Observations Y are MCAR (missing completely at random)
            mask <- runif(length(y)) < params$pr_missing_data
            y[mask] <- NA
        }

        return(list(field_state=fields[[field_id]]$state,
                    field_id=field_id,
                    y=y,
                    z=z,
                    time=time,
                    pixel_i=pixel_i,
                    pixel_j=pixel_j))
    })

    estimates <- get_hmm_and_minimum_distance_estimates_random_initialization(params=params0, panel=pixel_panel)

    return(list(pixel_panel=pixel_panel,
                estimates=estimates,
                params=params0,
                simulation_id=simulation_id))
}

## Simulation parameters
params0 <- get_params0()

## Overall average Pr[Y | S] used when include_z_in_simulation is false
params0$pr_y

## Values for Pr[Y | S, Z] used when include_z_in_simulation is true
params0$pr_y_given_cloudy <- rbind(c(0.82, 0.18),
                                   c(0.3, 0.7))
params0$pr_y_given_clear <- 2 * params0$pr_y - params0$pr_y_given_cloudy

## Sanity check
all(0.5 * params0$pr_y_given_clear + 0.5*params0$pr_y_given_cloudy == params0$pr_y)
all(params0$pr_y_given_clear <= 1)

## Observations are arranged in a square matrix of n_pixels_per_side^2 pixels in all simulations
n_pixels_per_side <- 100
adjacency <- adjacency.matrix(m=n_pixels_per_side, n=n_pixels_per_side)

## TODO Bump back up
n_simulations <- 50

for(pr_missing_data in c(0.0, 0.1)) {
    for(include_z_in_simulation in c(TRUE, FALSE)) {

        if(include_z_in_simulation) {
            ising_beta_values <- c(0.0, 2.0)
            ## Note: if Z is not constant over time, it is i.i.d. over time
            params0$z_constant_over_time <- FALSE
        } else {
            ## Note: the Ising beta parameter has no effect when we don't include Z in the simulation
            ising_beta_values <- c(0.0)
        }

        for(n_fields in c(10000, 100)) {
            for(ising_beta in ising_beta_values) {

                params0$include_z_in_simulation <- include_z_in_simulation
                params0$n_fields <- n_fields
                params0$ising_beta <- ising_beta
                params0$pr_missing_data <- pr_missing_data

                max_cores <- 12
                num_cores <- min(detectCores(), max_cores)
                cluster <- makeCluster(num_cores)  # Call stopCluster when done

                clusterExport(cluster, c("baum_welch",
                                         "eq_function_minimum_distance",
                                         "get_expectation_maximization_estimates",
                                         "get_hmm_and_minimum_distance_estimates_random_initialization",
                                         "get_min_distance_estimates",
                                         "get_random_initial_parameters",
                                         "get_transition_probs_from_M_S_joint",
                                         "is_diag_dominant",
                                         "objfn_minimum_distance",
                                         "run_single_simulation",
                                         "simulate_discrete_markov",
                                         "simulate_hmm",
                                         "simulate_ising",
                                         "valid_panel_element",
                                         "valid_parameters"))
            
                simulations <- parLapply(cluster,
                                         seq_len(n_simulations),
                                         run_single_simulation,
                                         params0=params0,
                                         adjacency=adjacency,
                                         n_pixels_per_side=n_pixels_per_side)
        
                stopCluster(cluster)

                z_correlation_description <- ""
                if(include_z_in_simulation) {
                    if(params0$z_constant_over_time) {
                        z_correlation_description <- "_constant_over_time"
                    } else {
                        z_correlation_description <- "_iid_over_time"
                    }
                }

                simulations_filename <- sprintf("spatial_corr_n_fields_%s_pr_missing_data_%s_include_z_%s%s_ising_beta_%s_%s_simulations.rds",
                                                params0$n_fields,
                                                pr_missing_data,
                                                include_z_in_simulation,
                                                z_correlation_description,
                                                params0$ising_beta,
                                                n_simulations)
                saveRDS(simulations, simulations_filename)

                ## TODO Include summary stats about Z in simulation summary?  Should be around 50-50 for high and low "clouds"
                ## TODO Include estimates of pr_y in simulation summary

                simulation_summaries <- lapply(simulations, function(simulation) {            
                    P_hat_naive <- lapply(simulation$estimates$M_Y_joint_hat, get_transition_probs_from_M_S_joint)
                    data.table(em_estimated_transition_1_2=sapply(simulation$estimates$em_params_hat_best_likelihood$P_list, function(x) return(x[1, 2])),
                               md_estimated_transition_1_2=sapply(simulation$estimates$min_dist_params_hat_best_objfn$P_list, function(x) return(x[1, 2])),
                               naive_estimated_transition_1_2=sapply(P_hat_naive, function(x) return(x[1, 2])),
                               true_transition_1_2=sapply(params0$P_list, function(x) return(x[1, 2])),
                               em_estimated_transition_2_1=sapply(simulation$estimates$em_params_hat_best_likelihood$P_list, function(x) return(x[2, 1])),
                               md_estimated_transition_2_1=sapply(simulation$estimates$min_dist_params_hat_best_objfn$P_list, function(x) return(x[2, 1])),
                               naive_estimated_transition_2_1=sapply(P_hat_naive, function(x) return(x[2, 1])),
                               true_transition_2_1=sapply(params0$P_list, function(x) return(x[2, 1])),
                               time=head(simulation$pixel_panel[[1]]$time, length(simulation$pixel_panel[[1]]$time) - 1),
                               simulation_id=simulation$simulation_id)
                })
    
                simulation_summary <- rbindlist(simulation_summaries)

                ## TODO Make sure this is working correctly when including both 1->2 and 2->1 transition probability estimates
                simulation_summary_melt  <- melt(simulation_summary, id.vars=c("simulation_id", "time", "true_transition_1_2", "true_transition_2_1"))
    
                simulation_summary_melt$time_label <- sprintf("time %s to %s", simulation_summary_melt$time, simulation_summary_melt$time+1)
    
                simulation_summary_melt[variable %like% "em_", algorithm := "ML"]
                simulation_summary_melt[variable %like% "md_", algorithm := "MD"]
                simulation_summary_melt[variable %like% "naive_", algorithm := "Frequency"]

                simulation_summary_melt[variable %like% "transition_1_2", state_label := "state 1 to 2"]
                simulation_summary_melt[variable %like% "transition_2_1", state_label := "state 2 to 1"]

                simulation_summary_melt[, true_transition := ifelse(state_label == "state 1 to 2", true_transition_1_2, true_transition_2_1)]
    
                simulation_summary_melt[, algorithm := factor(algorithm, levels=c("Frequency", "ML", "MD"))]

                field_width <- n_pixels_per_side / sqrt(n_fields)

                z_title_description <- "No Z"
                if(include_z_in_simulation) {
                    z_title_description <- sprintf("%s, Ising Beta = %s",
                                                   ifelse(params0$z_constant_over_time, "Z Constant Over Time", "Z I.I.D. Over Time"),
                                                   params0$ising_beta)
                }
                missing_data_description <- ifelse(pr_missing_data > 0, sprintf("Observations %s MCAR", pr_missing_data), "No Missing Data")
                title <- sprintf("%s Simulations, %s\n%s-by-%s Pixel Fields, %s",
                                 n_simulations,
                                 z_title_description,
                                 field_width,
                                 field_width,
                                 missing_data_description)

                p <- (ggplot(simulation_summary_melt, aes(y=value, x=algorithm, group=variable)) +
                      geom_boxplot() +
                      geom_hline(aes(yintercept=true_transition), linetype="dashed") +
                      ylab("transition probability") +
                      theme_bw() +
                      ggtitle(title) +
                      ylim(c(0, 0.65)) +
                      theme(plot.title = element_text(hjust = 0.5)) +
                      facet_grid(state_label ~ time_label))
                p
                filename <- sprintf("simulation_spatial_corr_%s_n_fields_%s_pr_missing_data_%s_include_z_%s%s_ising_beta_%s_%s_simulations.png",
                                    "estimated_transition_probabilities",
                                    params0$n_fields,
                                    pr_missing_data,
                                    include_z_in_simulation,
                                    z_correlation_description,
                                    params0$ising_beta,
                                    n_simulations)
                ggsave(filename, plot=p, width=6, height=4, units="in")

                ## These are plots showing a single simulation in detail
                dtable <- rbindlist(Map(data.frame, simulations[[1]]$pixel_panel))
                dtable[, time_label := sprintf("time %s", time)]

                dtable[, classification_error := field_state != y]
                dtable[, classification_error_label := ifelse(field_state != y, "Y ≠ S (misclassification)", "Y = S (correct classification)")]
                
                colors <- c("#dfc27d", "#018571")

                p <- (ggplot(dtable, aes(x=pixel_i, y=pixel_j, fill=factor(field_state))) +
                      facet_wrap(~ time_label) +
                      geom_raster() +
                      scale_fill_manual("true field state", values=colors) +
                      xlab("pixel coordinate (easting)") +
                      ylab("pixel coordinate (northing)"))

                filename <- sprintf("%s_n_fields_%s_pr_missing_data_%s_include_z_%s%s_ising_beta_%s_%s_simulations.png",
                                    "simulation_spatial_corr_true_field_state",
                                    params0$n_fields,
                                    pr_missing_data,
                                    include_z_in_simulation,
                                    z_correlation_description,
                                    params0$ising_beta,
                                    n_simulations)
                ggsave(filename, plot=p, width=6, height=4, units="in")

                title <- sprintf("Overall Misclassification Rate: %s", round(mean(dtable$classification_error, na.rm=TRUE), 3))
                p <- (ggplot(dtable, aes(x=pixel_i, y=pixel_j, fill=factor(classification_error_label))) +
                      facet_wrap(~ time_label) +
                      geom_raster() +
                      ggtitle(title) +
                      scale_fill_manual("classification outcome", values=c("grey", "red")) +
                      xlab("pixel coordinate (easting)") +
                      ylab("pixel coordinate (northing)"))

                filename <- sprintf("%s_n_fields_%s_pr_missing_data_%s_include_z_%s%s_ising_beta_%s_%s_simulations.png",
                                    "simulation_spatial_corr_classification_errors",
                                    params0$n_fields,
                                    pr_missing_data,
                                    include_z_in_simulation,
                                    z_correlation_description,
                                    params0$ising_beta,
                                    n_simulations)
                ggsave(filename, plot=p, width=6, height=4, units="in")

                p <- (ggplot(dtable, aes(x=pixel_i, y=pixel_j, fill=factor(z))) +
                      geom_raster() +
                      facet_wrap(~ time_label) +
                      scale_fill_manual("Z", values=c("grey", "black")) +
                      xlab("pixel coordinate (easting)") +
                      ylab("pixel coordinate (northing)"))

                filename <- sprintf("simulation_spatial_corr_z_n_fields_%s_pr_missing_data_%s_include_z_%s%s_ising_beta_%s_%s_simulations.png",
                                    params0$n_fields,
                                    pr_missing_data,
                                    include_z_in_simulation,
                                    z_correlation_description,
                                    params0$ising_beta,
                                    n_simulations)
                ggsave(filename, plot=p, width=6, height=4, units="in")

                p <- (ggplot(dtable, aes(x=pixel_i, y=pixel_j, fill=factor(y))) +
                      facet_wrap(~ time_label) +
                      geom_raster() +
                      scale_fill_manual("observed y", values=colors) +
                      xlab("pixel coordinate (easting)") +
                      ylab("pixel coordinate (northing)"))

                filename <- sprintf("simulation_spatial_corr_observed_y_n_fields_%s_pr_missing_data_%s_include_z_%s%s_ising_beta_%s_%s_simulations.png",
                                    params0$n_fields,
                                    pr_missing_data,
                                    include_z_in_simulation,
                                    z_correlation_description,
                                    params0$ising_beta,
                                    n_simulations)
                ggsave(filename, plot=p, width=6, height=4, units="in")
            }
        }
    }
}

## TODO Also save scatterplot showing frequency on x-axis, EM and MD on y-axis (or showing errors for all relative to true transition proba)
## TODO Also save plot of estimates of misclassification probabilities
