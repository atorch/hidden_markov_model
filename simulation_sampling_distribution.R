library(data.table)
library(ggplot2)
library(grid)
library(latex2exp)  # For ggplot2 xlab
library(parallel)
library(Rsolnp)

source("hmm_functions.R")
source("hmm_parameters.R")
source("ggplot_utils.R")

set_ggplot_theme()

set.seed(321321)

## True HMM parameters used in simulation
params0 <- get_params0()

get_estimates <- function(params0, n_panel_elements, n_random_starts) {

    panel <- replicate(n_panel_elements, simulate_hmm(params0), simplify=FALSE)
    estimates <- get_hmm_and_minimum_distance_estimates_random_initialization(params0, panel, n_random_starts)

    return(estimates)
}

cluster <- makeCluster(detectCores())
n_replications <- 100
n_random_starts <- 6
outfile_format <- "simulation_sampling_distribution_with_%s_random_initial_parameters_panel_size_%s_%s_replications.rds"
clusterExport(cluster, c("baum_welch",
                         "eq_function_minimum_distance",
                         "get_estimates",
                         "get_expectation_maximization_estimates",
                         "get_hmm_and_minimum_distance_estimates_random_initialization",
                         "get_min_distance_estimates",
                         "get_random_initial_parameters",
                         "get_transition_probs_from_M_S_joint",
                         "n_random_starts",
                         "objfn_minimum_distance",
                         "params0",
                         "simulate_discrete_markov",
                         "simulate_hmm",
                         "valid_panel_element",
                         "valid_parameters"))

panel_sizes <- c(100, 200, 500, 1000, 5000, 10000)  # Slow, panel size 5000 took 22 hours on my laptop (4 CPUs)
panel_sizes <- c(100, 200, 500, 1000)

for(panel_size in panel_sizes) {
    outfile <- sprintf(outfile_format, n_random_starts, panel_size, n_replications)
    if(file.exists(outfile)) {
        message(outfile, " already exists; won't re-run simulation")
    } else {
        message("running replications with random initialization, time is ", Sys.time())
        replications_random_initialization <- parLapply(cluster, rep(panel_size, n_replications), function(x) {
            get_estimates(params0=params0, n_panel_elements=x, n_random_starts=n_random_starts)  # Slow...
        })
        message("done, time is ", Sys.time())
        message("saving ", outfile)
        saveRDS(replications_random_initialization, outfile)
    }
}
stopCluster(cluster)

infiles <- sprintf(outfile_format, n_random_starts, panel_sizes, n_replications)

## This is a list of lists: for each set of parameter values, we have multiple replications of the simulation
all_simulations <- lapply(infiles, readRDS)
stopifnot(all(sapply(all_simulations, length) == n_replications))

dataframes <- lapply(all_simulations, function(replications) {
    df <- data.frame(panel_size=sapply(replications, function(replication) {
        return(replication$em_params_hat_best_likelihood$panel_size)
    }),
    em_params_hat_n_iterations=sapply(replications, function(replication) {
        return(replication$em_params_hat_best_likelihood$n_em_iterations)
    }),
    em_params_hat_pr_y_11=sapply(replications, function(replication) {
        return(replication$em_params_hat_best_likelihood$pr_y[1, 1])
    }),
    em_params_hat_pr_y_22=sapply(replications, function(replication) {
        return(replication$em_params_hat_best_likelihood$pr_y[2, 2])
    }),
    min_dist_params_hat_pr_y_11=sapply(replications, function(replication) {
        return(replication$min_dist_params_hat_best_objfn$pr_y[1, 1])
    }),
    min_dist_params_hat_pr_y_22=sapply(replications, function(replication) {
        return(replication$min_dist_params_hat_best_objfn$pr_y[2, 2])
    }))
    for(time_index in seq_along(params0$P_list)) {
        varname <- sprintf("em_params_hat_P%s_11", time_index)
        df[, varname] <- sapply(replications, function(replication) {
            return(replication$em_params_hat_best_likelihood$P_list[[time_index]][1, 1])
        })
        varname <- sprintf("min_dist_params_hat_P%s_11", time_index)
        df[, varname] <- sapply(replications, function(replication) {
            return(replication$min_dist_params_hat_best_objfn$P_list[[time_index]][1, 1])
        })
    }
    return(df)
})

df <- rbindlist(dataframes)
df$panel_size_label <- sprintf("%s%s", ifelse(df$panel_size %in% c(100, 1000), "sample size = ", ""), df$panel_size)
df$panel_size_label <- factor(df$panel_size_label,
                              levels=sprintf("%s%s", ifelse(sort(unique(df$panel_size)) %in% c(100, 1000), "sample size = ", ""), sort(unique(df$panel_size))))

estimators <- c("em", "min_dist")

## Transition probs
for (time_index in seq_along(params0$P_list)) {
    for(estimator in estimators) {
        variable <- sprintf("%s_params_hat_P%s_11", estimator, time_index)
        p <- (ggplot(df, aes_string(x=variable)) +
              geom_vline(xintercept=params0$P_list[[time_index]][1, 1], color="grey", linetype=2) +
              geom_histogram(binwidth=0.005, color="black", fill="white") +
              facet_wrap(~ panel_size_label, scale="free_x") +
              geom_rug(aes(x=true_value), data=data.frame(true_value=params0$P_list[[time_index]][1, 1])) +  # Correct parameter value
              ylab("") + xlab(""))
        outfile <- sprintf("simulation_methodology_%s_random_initialization_P_hat_%s_11_sampling_distribution.png", estimator, time_index)
        message("saving ", outfile)
        ggsave(outfile, p, width=12.5, height=8)
    }
}

## Misclassification probs
for(matrix_index in c("11", "22")) {
    for(estimator in estimators) {
        variable <- sprintf("%s_params_hat_pr_y_%s", estimator, matrix_index)
        true_value <- ifelse(matrix_index == "11", params0$pr_y[1, 1], params0$pr_y[2, 2])
        p <- (ggplot(df, aes_string(x=variable)) +
              geom_vline(xintercept=true_value, color="grey", linetype=2) +
              geom_histogram(binwidth=0.005, color="black", fill="white") +
              facet_wrap(~ panel_size_label, scale="free_x") +
              geom_rug(aes(x=true_value), data=data.frame(true_value=true_value)) +  # Correct parameter value
              ylab("") + xlab(""))
        outfile <- sprintf("simulation_methodology_%s_random_initialization_pr_y_%s_sampling_distribution.png", estimator, matrix_index)
        message("saving ", outfile)
        ggsave(outfile, p, width=12.5, height=8)
    }
}

df_melted <- melt(df, id.vars=c("panel_size", "panel_size_label"))
p <- (ggplot(df_melted, aes(x=value)) +
      geom_histogram(binwidth=0.01, fill="white", color="black") +
      facet_wrap(~ variable, scale="free_x"))
p  # Looks correct, compare to params0 -- careful, lots of estimates hit max_iter

## Do min dist estimates vary with initial parameters?
P_hat_range_min_dist <- c(lapply(all_simulations, function(list_of_replications) {
    sapply(list_of_replications, function(replication) {
        matrix_of_P_hat <- sapply(replication$min_dist_params_hat, function(x) {
            return(c(x$P_list, recursive=TRUE))
        })
        P_hat_range <- max(apply(matrix_of_P_hat, 1, max) - apply(matrix_of_P_hat, 1, min))
        return(P_hat_range)
    })
}), recursive=TRUE)
sum(P_hat_range_min_dist > 0.01)  # 60 simulations out of length(P_hat_range_min_dist)=1200 where P_hat differs by more than 0.01 across initializations
sum(P_hat_range_min_dist > 0.10)  # 37 simulations out of length(P_hat_range_min_dist)=1200 where P_hat differs by more than 0.10 across initializations
P_hat_range_em <- c(lapply(all_simulations, function(list_of_replications) {
    sapply(list_of_replications, function(replication) {
        matrix_of_P_hat <- sapply(replication$em_params_hat_list, function(x) {
            return(c(x$P_list, recursive=TRUE))
        })
        P_hat_range <- max(apply(matrix_of_P_hat, 1, max) - apply(matrix_of_P_hat, 1, min))
        return(P_hat_range)
    })
}), recursive=TRUE)
sum(P_hat_range_em > 0.01)  # 1183 simulations out of length(P_hat_range_em)=1200
sum(P_hat_range_em > 0.10)  # 394 simulations out of length(P_hat_range_em)=1200

## RMSE of each estimator as fn of sample size
dtable <- data.table(df)
setkey(dtable, panel_size)
with(subset(dtable, panel_size <= 5000), t.test(x=(em_params_hat_pr_y_11 - params0$pr_y[1, 1])^2,
                                                y=(min_dist_params_hat_pr_y_11 - params0$pr_y[1, 1])^2))  # P-value 0.02
with(subset(dtable, panel_size <= 5000), t.test(x=(em_params_hat_P1_11 - params0$P_list[[1]][1, 1])^2,
                                                y=(min_dist_params_hat_P1_11 - params0$P_list[[1]][1, 1])^2))  # P-value 0.11
with(subset(dtable, panel_size <= 5000), t.test(x=(em_params_hat_P2_11 - params0$P_list[[2]][1, 1])^2,
                                                y=(min_dist_params_hat_P2_11 - params0$P_list[[2]][1, 1])^2))  # P-value 0.09
with(subset(dtable, panel_size <= 5000), t.test(x=(em_params_hat_P3_11 - params0$P_list[[3]][1, 1])^2,
                                                y=(min_dist_params_hat_P3_11 - params0$P_list[[3]][1, 1])^2))  # P-value 0.004
dtable_rmses <- dtable[, list(rmse_em_pr_y_11=sqrt(mean((em_params_hat_pr_y_11 - params0$pr_y[1, 1])^2)),
                              rmse_min_dist_pr_y_11=sqrt(mean((min_dist_params_hat_pr_y_11 - params0$pr_y[1, 1])^2, na.rm=TRUE)),
                              rmse_em_P1_11=sqrt(mean((em_params_hat_P1_11 - params0$P_list[[1]][1, 1])^2)),
                              rmse_min_dist_P1_11=sqrt(mean((min_dist_params_hat_P1_11 - params0$P_list[[1]][1, 1])^2, na.rm=TRUE)),
                              rmse_em_P2_11=sqrt(mean((em_params_hat_P2_11 - params0$P_list[[2]][1, 1])^2)),
                              rmse_min_dist_P2_11=sqrt(mean((min_dist_params_hat_P2_11 - params0$P_list[[2]][1, 1])^2, na.rm=TRUE)),
                              rmse_em_P3_11=sqrt(mean((em_params_hat_P3_11 - params0$P_list[[3]][1, 1])^2)),
                              rmse_min_dist_P3_11=sqrt(mean((min_dist_params_hat_P3_11 - params0$P_list[[3]][1, 1])^2, na.rm=TRUE))),
                       by="panel_size"]
dtable_melted <- melt(dtable_rmses, id.vars="panel_size")
dtable_melted[, estimator := ifelse(grepl("_em_", variable), "EM", "MinDist")]
dtable_melted[, parameter := ifelse(grepl("pr_y_11", variable), "first entry of misclassification matrix",
                                    sprintf("first entry of %s", gsub("rmse_min_dist_|rmse_em_|_11", "", variable)))]

p <- (ggplot(dtable_melted, aes(x=panel_size, y=value, color=estimator)) +
      geom_point() + geom_line(size=1.2) +
      scale_color_manual("estimator", values=c("black", "grey")) +
      scale_x_continuous("panel size", breaks=c(100, 1000, 5000, 10000)) +
      scale_y_continuous("root mean squared error") +
      geom_hline(yintercept=0, linetype=2, color="grey") +
      facet_wrap(~ parameter, scales="free_y"))
p
outfile <- "simulation_methodology_random_initialization_rmse_em_and_min_dist.png"
message("saving ", outfile)
ggsave(outfile, p, width=11, height=8)

message("done with simulation, time is ", Sys.time())
