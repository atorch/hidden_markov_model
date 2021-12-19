library(ggplot2)
library(rgdal)
library(rgeos)
library(sp)

fileDir <- "."
pattern <- "estimates_window_[0-9]*00[01]_[0-9]*00[01]_width_1000_class_frequency_cutoff_0.005_subsample_0.01_combined_classes_grassland_as_forest_combine_other_non_forest_use_md_as_initial_values_for_em.rds"
estimate_filenames <- list.files(path=fileDir,
                                 pattern=pattern,
                                 full.names=TRUE)

crs_longlat <- CRS("+proj=longlat")

brazil_states <- readOGR(dsn="./state_boundaries/", layer="BRUFE250GC_SIR")
brazil_state_names <- brazil_states$NM_ESTADO
brazil_state_polygons <- as(brazil_states, "SpatialPolygons")

brazil_municipalities <- readOGR(dsn="./munic_boundaries/", layer="BRMUE250GC_SIR")
brazil_municipality_names <- brazil_municipalities$NM_MUNICIP
brazil_municipality_polygons <- as(brazil_municipalities, "SpatialPolygons")

# TODO Use spTransform here?  These are originally "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"
proj4string(brazil_state_polygons) <- crs_longlat
proj4string(brazil_municipality_polygons) <- crs_longlat

forest_class <- 3
agriculture_and_pasture_class <- 21

## Compute the marginal distribution over hidden states
calculate_mu_t <- function(time_index, params) {
    if(time_index == 1) {
        mu_t <- params$mu  # Equals initial distribution when t=1
    } else {
        mu_t <- params$mu %*% Reduce("%*%", params$P_list[seq_len(time_index - 1)])
    }
    stopifnot(isTRUE(all.equal(sum(mu_t), 1)))  # Valid probability distribution, careful comparing floats
    return(mu_t)
}

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

        df <- data.frame(deforestation_rate_ml=1 - pr_remain_forest_ml,
                         deforestation_rate_md=1 - pr_remain_forest_md,
                         deforestation_rate_freq=1 - pr_remain_forest_freq,
                         time_index=seq_along(pr_remain_forest_ml),
                         window_row=estimates$options$row,
                         window_col=estimates$options$col)

        df$hidden_state_mapbiomas_classes <- paste(estimates$mapbiomas_classes_to_keep, collapse=",")

        df$pr_y_diagonal_sum_ml <- sum(diag(estimates$em_params_hat_best_likelihood$pr_y))
        df$pr_y_diagonal_sum_md <- sum(diag(estimates$min_dist_params_hat_best_objfn$pr_y))
        
        df$pr_y_diagonal_forest_ml <- estimates$em_params_hat_best_likelihood$pr_y[forest_index, forest_index]
        df$pr_y_diagonal_forest_md <- estimates$min_dist_params_hat_best_objfn$pr_y[forest_index, forest_index]

        df$fraction_missing_in_all_years <- estimates$fraction_missing_in_all_years

        mu_t_ml <- lapply(df$time_index, calculate_mu_t, params=estimates$em_params_hat_best_likelihood)
        df$fraction_forest_ml <- sapply(mu_t_ml, function(mu) mu[forest_index])

        df$reforestation_rate_ml <- sapply(df$time_index, function(t) {
            weights <- mu_t_ml[[t]][-forest_index] / sum(mu_t_ml[[t]][-forest_index])
            return(sum(weights * estimates$em_params_hat_best_likelihood$P_list[[t]][-forest_index, forest_index]))
        })

        ## Sanity check: fraction_forest_ml should be consistent with the deforestation and reforestation rates
        stopifnot(isTRUE(all.equal(df$fraction_forest_ml[2],
        (1 - df$fraction_forest_ml[1]) * df$reforestation_rate_ml[1] + df$fraction_forest_ml[1] * (1 - df$deforestation_rate_ml[1]))))

        mu_t_md <- lapply(df$time_index, calculate_mu_t, params=estimates$min_dist_params_hat_best_objfn)
        df$fraction_forest_md <- sapply(mu_t_md, function(mu) mu[forest_index])

        df$reforestation_rate_md <- sapply(df$time_index, function(t) {
            weights <- mu_t_md[[t]][-forest_index] / sum(mu_t_md[[t]][-forest_index])
            return(sum(weights * estimates$em_params_hat_best_likelihood$P_list[[t]][-forest_index, forest_index]))
        })

        if(agriculture_and_pasture_class %in% estimates$mapbiomas_classes_to_keep) {
            agriculture_and_pasture_index <- which(estimates$mapbiomas_classes_to_keep == agriculture_and_pasture_class)
            df$pr_agriculture_and_pasture_to_forest_ml <- sapply(estimates$em_params_hat_best_likelihood$P_list,
                                                                 function(P) P[agriculture_and_pasture_index, forest_index])
            df$pr_agriculture_and_pasture_to_forest_md <- sapply(estimates$min_dist_params_hat_best_objfn$P_list,
                                                                 function(P) P[agriculture_and_pasture_index, forest_index])
            df$pr_agriculture_and_pasture_to_forest_freq <- sapply(estimates$P_hat_frequency,
                                                                   function(P) P[agriculture_and_pasture_index, forest_index])
            df$fraction_agriculture_and_pasture_ml <- sapply(mu_t_ml, function(mu) mu[agriculture_and_pasture_index])
            df$fraction_agriculture_and_pasture_md <- sapply(mu_t_md, function(mu) mu[agriculture_and_pasture_index])
        } else {
            ## These are undefined in windows without agriculture and pasture
            df$pr_agriculture_and_pasture_to_forest_ml <- NA
            df$pr_agriculture_and_pasture_to_forest_md <- NA
            df$pr_agriculture_and_pasture_to_forest_freq <- NA
            df$fraction_agriculture_and_pasture_ml <- NA
            df$fraction_agriculture_and_pasture_md <- NA
        }

        ## Fraction of forest classifications _after_ combining classes, all years combined
        df$forest_frequency_all_years <- estimates$class_frequencies[names(estimates$class_frequencies) == forest_class]

        window_bbox <- bbox(t(array(estimates$window_bbox)))
        window_polygon <- as(raster::extent(window_bbox), "SpatialPolygons")
        window_centroid <- coordinates(window_polygon)

        proj4string(window_polygon) <- crs_longlat

        intersected_states <- brazil_state_names[over(window_polygon, brazil_state_polygons)]
        intersected_municipalities <- brazil_municipality_names[over(window_polygon, brazil_municipality_polygons)]

        df$states <- paste(intersected_states, collapse=",")
        df$municipalities <- paste(intersected_municipalities, collapse=",")

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

df$window_id <- sprintf("%s_%s", df$window_row, df$window_col)
df$year = df$time_index + 1984

filename <- sprintf("estimated_deforestation_rates_%s_use_md_as_initial_values_for_em.csv", format(Sys.time(), "%Y_%m_%d"))
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
