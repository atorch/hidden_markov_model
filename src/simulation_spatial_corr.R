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

    ## TODO Might be cleaner to have z in {1, 2} instead of in {-1, 1},
    ## so that it can be used as an array index
    ## TODO Speed this up, and then increase n_iter

    df$z <- simulate_ising(n_pixels=nrow(df), adjacency=adjacency, beta=params0$ising_beta, n_iter=50)

    n_fields_per_side <- sqrt(n_fields)

    ## TODO Could have fields of different sizes by changing the edges here
    ## TODO Could speed up simulation by doing this once outside this fn and passing it in
    field_edges <- expand.grid(pixel_i=seq(1, n_pixels_per_side, n_pixels_per_side / n_fields_per_side),
                               pixel_j=seq(1, n_pixels_per_side, n_pixels_per_side / n_fields_per_side))
    stopifnot(nrow(field_edges) == n_fields)

    for(field_id in seq_len(n_fields)) {
        cutoff_i <- field_edges[field_id, ]$pixel_i
        cutoff_j <- field_edges[field_id, ]$pixel_j
        df[df$pixel_i >= cutoff_i & df$pixel_j >= cutoff_j, ]$field_id <- field_id
    }

    pixel_panel <- lapply(seq_len(nrow(df)), function(row_index, params=params0) {
        field_id <- df[row_index, ]$field_id
        pixel_i <- df[row_index, ]$pixel_i
        pixel_j <- df[row_index, ]$pixel_j
        z <- df[row_index, ]$z
        field_state <- fields[[field_id]]$state
        y <- vapply(field_state, function(s) {
            if(z == 1) {
                sample(seq_len(ncol(params$pr_y)), size=1, prob=params$pr_y_cloudy[s, ])
            } else {
                sample(seq_len(ncol(params$pr_y)), size=1, prob=params$pr_y_given_clear[s, ])
            }
        }, FUN.VALUE=1)
        time <- seq_along(y)
        return(list(field_state=field_state,
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

## Overall average Pr[Y | S]
params0$pr_y
params0$pr_y_given_cloudy <- rbind(c(0.82, 0.18),
                                   c(0.3, 0.7))
params0$pr_y_given_clear <- 2 * params0$pr_y - params0$pr_y_given_cloudy

## Sanity check
all(0.5 * params0$pr_y_given_clear + 0.5*params0$pr_y_given_cloudy == params0$pr_y)
all(params0$pr_y_given_clear <= 1)

n_pixels_per_side <- 100
adjacency <- adjacency.matrix(m=n_pixels_per_side, n=n_pixels_per_side)

n_simulations <- 100

for(n_fields in c(10000, 400, 100)) {
    for(ising_beta in c(0.0, 0.5, 1.0, 1.5, 2.0)) {

        params0$n_fields <- n_fields
        params0$ising_beta <- ising_beta        

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

        simulations_filename <- sprintf("spatial_corr_n_fields_%s_ising_beta_%s_%s_simulations.rds",
                                        params0$n_fields,
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

        p <- (ggplot(simulation_summary, aes(y=em_estimated_transition_1_2, x=naive_estimated_transition_1_2)) +
              geom_point() +
              theme_bw())
        p
    
        p <- (ggplot(simulation_summary, aes(y=md_estimated_transition_1_2, x=naive_estimated_transition_1_2)) +
              geom_point() +
              theme_bw())
        p

        ## TODO Make sure this is working correctly when including both 1->2 and 2->1 transition probability estimates
        simulation_summary_melt  <- melt(simulation_summary, id.vars=c("simulation_id", "time", "true_transition_1_2", "true_transition_2_1"))
    
        simulation_summary_melt$time_label <- sprintf("time %s to %s", simulation_summary_melt$time, simulation_summary_melt$time+1)
    
        simulation_summary_melt[variable %like% "em_", algorithm := "EM"]
        simulation_summary_melt[variable %like% "md_", algorithm := "MD"]
        simulation_summary_melt[variable %like% "naive_", algorithm := "Frequency"]

        simulation_summary_melt[variable %like% "transition_1_2", state_label := "state 1 to 2"]
        simulation_summary_melt[variable %like% "transition_2_1", state_label := "state 2 to 1"]

        simulation_summary_melt[, true_transition := ifelse(state_label == "state 1 to 2", true_transition_1_2, true_transition_2_1)]
    
        simulation_summary_melt[, algorithm := factor(algorithm, levels=c("Frequency", "EM", "MD"))]

        field_width <- n_pixels_per_side / sqrt(n_fields)

        title <- sprintf("%s Simulations, Ising Beta = %s, %s-by-%s Pixel Fields", n_simulations, params0$ising_beta, field_width, field_width)
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
        filename <- sprintf("simulation_spatial_corr_estimated_transition_probabilities_n_fields_%s_ising_beta_%s_%s_simulations.png",
                            params0$n_fields,
                            params0$ising_beta,
                            n_simulations)
        ggsave(filename, plot=p, width=6, height=4, units="in")

        ## These are plots showing a single simulation in detail
        dtable <- rbindlist(Map(data.frame, simulations[[1]]$pixel_panel))
        dtable[, time_label := sprintf("time %s", time)]
        dtable[, classification_error := field_state != y]

        colors <- c("#dfc27d", "#018571")

        p <- (ggplot(dtable, aes(x=pixel_i, y=pixel_j, fill=factor(field_state))) +
              facet_wrap(~ time_label) +
              geom_raster() +
              scale_fill_manual("true field state", values=colors) +
              xlab("pixel coordinate (easting)") +
              ylab("pixel coordinate (northing)"))

        filename <- sprintf("simulation_spatial_corr_true_field_state_n_fields_%s_ising_beta_%s_%s_simulations.png",
                            params0$n_fields,
                            params0$ising_beta,
                            n_simulations)
        ggsave(filename, plot=p, width=6, height=4, units="in")

        p <- (ggplot(dtable, aes(x=pixel_i, y=pixel_j, fill=factor(classification_error))) +
              facet_wrap(~ time_label) +
              geom_raster() +
              scale_fill_manual("classification error", values=c("grey", "red")) +
              xlab("pixel coordinate (easting)") +
              ylab("pixel coordinate (northing)"))
        filename <- sprintf("simulation_spatial_corr_classification_errors_n_fields_%s_ising_beta_%s_%s_simulations.png",
                            params0$n_fields,
                            params0$ising_beta,
                            n_simulations)
        ggsave(filename, plot=p, width=6, height=4, units="in")

        ## TODO Save plot of Z -- static or changing over time?
        p <- (ggplot(dtable, aes(x=pixel_i, y=pixel_j, fill=factor(z))) +
              geom_raster() +
              facet_wrap(~ time_label) +
              scale_fill_manual("z (aka 'clouds'\nor 'haze')", values=c("grey", "black")) +
              xlab("pixel coordinate (easting)") +
              ylab("pixel coordinate (northing)"))
        filename <- sprintf("simulation_spatial_corr_z_static_clouds_haze_n_fields_%s_ising_beta_%s_%s_simulations.png",
                            params0$n_fields,
                            params0$ising_beta,
                            n_simulations)
        ggsave(filename, plot=p, width=6, height=4, units="in")

        p <- (ggplot(dtable, aes(x=pixel_i, y=pixel_j, fill=factor(y))) +
              facet_wrap(~ time_label) +
              geom_raster() +
              scale_fill_manual("observed y", values=colors) +
              xlab("pixel coordinate (easting)") +
              ylab("pixel coordinate (northing)"))
        filename <- sprintf("simulation_spatial_corr_observed_y_n_fields_%s_ising_beta_%s_%s_simulations.png",
                            params0$n_fields,
                            params0$ising_beta,
                            n_simulations)
        ggsave(filename, plot=p, width=6, height=4, units="in")
    }
}

## TODO Also save scatterplot showing frequency on x-axis, EM and MD on y-axis (or showing errors for all relative to true transition proba)
## TODO Also save plot of estimates of misclassification probabilities

