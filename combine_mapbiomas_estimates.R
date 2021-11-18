library(ggplot2)

estimate_filenames <- list.files(pattern="estimates_window_[0-9]*_[0-9]*_width_500_class_frequency_cutoff_0.005_subsample_0.1_combined_classes_grassland_as_forest.rds")

forest_class <- 3

estimate_dfs <- list()
for(filename in estimate_filenames) {
    message("Loading ", filename)
    estimates <- readRDS(filename)
    if(forest_class %in% estimates$mapbiomas_classes_to_keep) {
        forest_index <- which(estimates$mapbiomas_classes_to_keep == forest_class)
        pr_remain_forest_ml <- sapply(estimates$em_params_hat_best_likelihood$P_list, function(P) P[forest_index, forest_index])
        pr_remain_forest_md <- sapply(estimates$min_dist_params_hat_best_objfn$P_list, function(P) P[forest_index, forest_index])
        pr_remain_forest_freq <- sapply(estimates$P_hat_frequency, function(P) P[forest_index, forest_index])
        df <- data.frame(deforestation_rate_ml=1 - pr_remain_forest_ml,
                         deforestation_rate_md=1 - pr_remain_forest_md,
                         deforestation_rate_freq=1 - pr_remain_forest_freq,
                         time_index=seq_along(pr_remain_forest_ml),
                         window_row=estimates$options$row,
                         window_col=estimates$options$col)

        ## TODO Get window lat lon from window_bbox, if necessary
        df$pr_y_diag_dominant_ml <- all(diag(estimates$em_params_hat_best_likelihood$pr_y) > 0.51)
        df$pr_y_diag_dominant_md <- all(diag(estimates$min_dist_params_hat_best_objfn$pr_y) > 0.51)  # Might be exactly 0.50 (hits constraint)
        df$n_mapbiomas_classes <- length(estimates$mapbiomas_classes_to_keep)

        estimate_dfs[[length(estimate_dfs) + 1]] <- df
    } else {
        message(filename, " doesn't contain class ", forest_class, ", skipping")
    }
}

df <- do.call(rbind, estimate_dfs)
filename <- sprintf("estimated_deforestation_rates_%s.csv", format(Sys.time(), "%Y_%m_%d"))
message("Writing ", filename, ", dataframe dim is ", nrow(df), " by ", ncol(df))
write.csv(df, filename, row.names=FALSE)

p <- (ggplot(df, aes(x=deforestation_rate_freq, y=deforestation_rate_ml)) +
      geom_abline(slope=1, lty=2, alpha=0.5) +
      geom_point(alpha=0.5) +
      theme_bw())
filename <- "estimated_deforestation_rates.png"
ggsave(p, filename=filename, width=6, height=4, units="in")

p <- (ggplot(subset(df, pr_y_diag_dominant_ml), aes(x=deforestation_rate_freq, y=deforestation_rate_ml)) +
      geom_abline(slope=1, lty=2, alpha=0.5) +
      geom_point(alpha=0.25) +
      ggtitle("Windows where estimated Pr[ Y | S ] is diag dominant") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)))
filename <- "estimated_deforestation_rates_diag_dominant.png"
ggsave(p, filename=filename, width=6, height=4, units="in")