library(ggplot2)
library(sp)

fileDir <- "."
pattern <- "estimates_window_[0-9]*00[01]_[0-9]*00[01]_width_1000_class_frequency_cutoff_0.005_subsample_0.01_combined_classes_grassland_as_forest_combine_other_non_forest_skip_ml_if_md_is_diag_dominant.rds"
estimate_filenames <- list.files(path=fileDir,
                                 pattern=pattern,
                                 full.names=TRUE)

forest_class <- 3

message("Found ", length(estimate_filenames), " .rds files matching pattern")

estimate_dfs <- list()
for(filename in estimate_filenames) {
    message("Loading ", filename)
    estimates <- readRDS(filename)
    if(forest_class %in% estimates$mapbiomas_classes_to_keep) {
        forest_index <- which(estimates$mapbiomas_classes_to_keep == forest_class)
        pr_remain_forest_md <- sapply(estimates$min_dist_params_hat_best_objfn$P_list, function(P) P[forest_index, forest_index])
        if("em_params_hat_best_likelihood" %in% names(estimates)) {
            pr_remain_forest_ml <- sapply(estimates$em_params_hat_best_likelihood$P_list, function(P) P[forest_index, forest_index])
        } else {
            ## If using --skip_ml_if_md_is_diag_dominant, this will happen in windows where MD's pr_y is diagonally dominant
            pr_remain_forest_ml <- rep(NA, length(pr_remain_forest_md))
        }        
        pr_remain_forest_freq <- sapply(estimates$P_hat_frequency, function(P) P[forest_index, forest_index])

        ## TODO We could also add a column for the probability of transitioning from forest to {crops, pasture}
        df <- data.frame(deforestation_rate_ml=1 - pr_remain_forest_ml,
                         deforestation_rate_md=1 - pr_remain_forest_md,
                         deforestation_rate_freq=1 - pr_remain_forest_freq,
                         time_index=seq_along(pr_remain_forest_ml),
                         window_row=estimates$options$row,
                         window_col=estimates$options$col)

        window_bbox <- bbox(t(array(estimates$window_bbox)))
        window_polygon <- as(raster::extent(window_bbox), "SpatialPolygons")
        window_centroid <- coordinates(window_polygon)

        ## TODO We could use either window_{lat, lon} or window_polygon to intersect with states/municipalities
        df$window_lat <- window_centroid[2]
        df$window_lon <- window_centroid[1]

        df$pr_y_diag_dominant_ml <- all(diag(estimates$em_params_hat_best_likelihood$pr_y) > 0.51)
        df$pr_y_diag_dominant_md <- all(diag(estimates$min_dist_params_hat_best_objfn$pr_y) > 0.51)  # Might be exactly 0.50 (hits constraint)
        df$n_mapbiomas_classes <- length(estimates$mapbiomas_classes_to_keep)

        estimate_dfs[[length(estimate_dfs) + 1]] <- df
    } else {
        message(filename, " doesn't contain class ", forest_class, ", skipping")
    }
}

df <- do.call(rbind, estimate_dfs)

message("Fraction of windows with diag dominant Pr[ Y | S ] for EM/ML:")
print(mean(df$pr_y_diag_dominant_ml))

message("Fraction of windows with diag dominant Pr[ Y | S ] for MD:")
print(mean(df$pr_y_diag_dominant_md))

message("Fraction of windows with diag dominant Pr[ Y | S ] for MD or EM/ML:")
print(mean(df$pr_y_diag_dominant_md | df$pr_y_diag_dominant_ml))

print(table(df$n_mapbiomas_classes))
print(table(df$n_mapbiomas_classes, df$pr_y_diag_dominant_md))

message("Summary stats for EM/ML deforestation rates:")
print(summary(df$deforestation_rate_ml))
message("Summary stats for MD deforestation rates:")
print(summary(df$deforestation_rate_md))

filename <- sprintf("estimated_deforestation_rates_%s.csv", format(Sys.time(), "%Y_%m_%d"))
message("Writing ", filename, ", dataframe dim is ", nrow(df), " by ", ncol(df))
write.csv(df, filename, row.names=FALSE)

title <- sprintf("Correlation = %s", round(cor(df$deforestation_rate_ml, df$deforestation_rate_md), 3))
p <- (ggplot(df, aes(x=deforestation_rate_ml, y=deforestation_rate_md)) +      
      geom_point(alpha=0.15) +
      geom_smooth(method="lm", formula=y ~ x) +
      geom_abline(slope=1, lty=2, alpha=0.5) +      
      ggtitle(title) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)))
filename <- "estimated_deforestation_rates_md_versus_ml.png"
ggsave(p, filename=filename, width=6, height=4, units="in")

p <- (ggplot(df, aes(x=deforestation_rate_freq, y=deforestation_rate_ml)) +
      geom_abline(slope=1, lty=2, alpha=0.5) +
      geom_point(alpha=0.15) +
      theme_bw())
filename <- "estimated_deforestation_rates_ml.png"
ggsave(p, filename=filename, width=6, height=4, units="in")

p <- (ggplot(subset(df, pr_y_diag_dominant_ml), aes(x=deforestation_rate_freq, y=deforestation_rate_ml)) +
      geom_abline(slope=1, lty=2, alpha=0.5) +
      geom_point(alpha=0.15) +
      ggtitle("Windows where estimated Pr[ Y | S ] is diag dominant") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)))
filename <- "estimated_deforestation_rates_ml_diag_dominant.png"
ggsave(p, filename=filename, width=6, height=4, units="in")

p <- (ggplot(df, aes(x=deforestation_rate_freq, y=deforestation_rate_md)) +
      geom_abline(slope=1, lty=2, alpha=0.5) +
      geom_point(alpha=0.15) +
      theme_bw())
filename <- "estimated_deforestation_rates_md.png"
ggsave(p, filename=filename, width=6, height=4, units="in")

p <- (ggplot(subset(df, pr_y_diag_dominant_md), aes(x=deforestation_rate_freq, y=deforestation_rate_md)) +
      geom_abline(slope=1, lty=2, alpha=0.5) +
      geom_point(alpha=0.15) +
      ggtitle("Windows where estimated Pr[ Y | S ] is diag dominant") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)))
filename <- "estimated_deforestation_rates_md_diag_dominant.png"
ggsave(p, filename=filename, width=6, height=4, units="in")
