library(Rsolnp)
library(data.table)
library(ggplot2)
library(optparse)
library(raster)

source("hmm_functions.R")

opt_list <- list(make_option("--mapbiomas_raster_path", default="./HMM_MapBiomas_v2/mapbiomas.vrt"),
                 make_option("--row", default=50000, type="integer"),
                 make_option("--col", default=51000, type="integer"),
                 make_option("--width_in_pixels", default=1000, type="integer"),
                 make_option("--subsample", default=0.1, type="double"),
                 make_option("--class_frequency_cutoff", default=0.005, type="double"),
                 make_option("--n_random_starts_em", default=2, type="integer"),
                 make_option("--n_random_starts_md", default=1, type="integer"),
                 make_option("--grassland_as_forest", default=FALSE, action="store_true"),
                 make_option("--combine_other_non_forest", default=FALSE, action="store_true"),
                 make_option("--skip_ml_if_md_is_diag_dominant", default=FALSE, action="store_true"),
                 make_option("--use_md_as_initial_values_for_em", default=FALSE, action="store_true"))
opt <- parse_args(OptionParser(option_list=opt_list))
message("command line options: ", paste(sprintf("%s=%s", names(opt), opt), collapse=", "))

mapbiomas <- stack(opt$mapbiomas_raster_path)
nlayers(mapbiomas)

window <- getValuesBlock(mapbiomas,
                         row=opt$row,
                         col=opt$col,
                         nrows=opt$width_in_pixels,
                         ncols=opt$width_in_pixels)

dim(window)
window_extent <- extent(mapbiomas, opt$row, opt$row + opt$width_in_pixels, opt$col, opt$col + opt$width_in_pixels)

window_raster<- raster(window_extent, crs=crs(mapbiomas), nrows=opt$width_in_pixels, ncols=opt$width_in_pixels)
for(time_index in c(1, 8)) {
    values(window_raster) <- window[, time_index]
    filename <- sprintf("./atlantic_forest_output/raster_window_%s_%s_width_%s_band_%s.tif", opt$row, opt$col, opt$width_in_pixels, time_index)

    ## These .tifs aren't used anywhere in the code, but it can be helpful to inspect these rasters in qgis
    message("Writing ", filename)
    writeRaster(window_raster, filename, overwrite=TRUE)
}

class_frequencies_before_combining <- round(table(window) / (nrow(window) * ncol(window)), 4)

pr_missing <- mean(is.na(window))
pr_water_or_sand <- mean(window %in% c(22, 33))
if(pr_missing > 0.9 || pr_water_or_sand > 0.5) {
    message("Window ", opt$row, " ", opt$col, " is missing at rate ",
            pr_missing, ", ",
            pr_water_or_sand, " water or sand (averaging over all bands), ",
            "skipping estimation")
    quit()
}

n_years <- ncol(window)
for(time_index in seq_len(n_years)) {
    pr_missing <- mean(is.na(window[, time_index]))
    if(pr_missing > 0.9) {
        message("Window ", opt$row, " ", opt$col, " is missing at rate ",
                pr_missing, " at time index (i.e. band) ", time_index,
                ", skipping estimation")
        quit()
    }
}

fraction_missing_in_all_years <- mean(rowMeans(is.na(window)) == 1.0)
count_missing_in_all_years <- sum(rowMeans(is.na(window)) == 1.0)
message("Fraction of pixels missing in 100% of years in the original data: ", fraction_missing_in_all_years)

## When constructing our panel (for estimation), we will only consider pixels that contain at least one non-missing observation
## in the original data. This will remove pixels in the ocean and pixels outside of the Atlantic forest region
valid_pixel_index <- rowMeans(is.na(window)) < 1.0

## Combine classes
## Class 12 (grassland) is optionally combined with class 3 (forest)
if(opt$grassland_as_forest) window[window %in% 12] <- 3

## Combine classes
## Classes 4 (savanna formation) and 9 (forest plantation) are combined with class 3 (forest)
window[window %in% c(4, 9)] <- 3

## Combine classes
## Class 11 (wetlands), class 22 (sand), and class 29 (rocky outcrop) are combined with class 33 (rivers and lakes)
window[window %in% c(11, 22, 29)] <- 33

## Combine classes
## Class 13 (other non-forest) is combined with class 33 (already a combination of wetlands, sand, rivers and lakes)
if(opt$combine_other_non_forest) window[window %in% 13] <- 33

## See https://mapbiomas-br-site.s3.amazonaws.com/downloads/Colecction%206/Cod_Class_legenda_Col6_MapBiomas_BR.pdf
unique_mapbiomas_classes <- sort(unique(c(window, recursive=TRUE)))

rare_mapbiomas_classes <- vector("numeric")
for(class in unique_mapbiomas_classes) {
    if(mean(window == class, na.rm=TRUE) < opt$class_frequency_cutoff) {
        rare_mapbiomas_classes <- c(rare_mapbiomas_classes, class)
    }
}

## We are going to recode rare classes as NA
## This is effectively assuming that all observations of rare classes must be misclassifications
## In most windows we will keep classes 3 and 21 (forest and mosaic of pasture + agriculture)
mapbiomas_classes_to_keep <- unique_mapbiomas_classes[!unique_mapbiomas_classes %in% rare_mapbiomas_classes]
message("Keeping the following classes:")
print(mapbiomas_classes_to_keep)

table(window)
class_frequencies <- round(table(window) / (nrow(window) * ncol(window)), 4)

## Careful, there can be missing values (even before we recode rare classes as NA)!
message("Missing value counts in the original data (fraction and sum):")
mean(is.na(c(window, recursive=TRUE)))
sum(is.na(c(window, recursive=TRUE)))

## Note: the code expects observations to be in the set {1, 2, 3, ..., |Y|},
## So we need to recode sets of classes like {3, 21} to {1, 2} for example
window_recoded <- window

for(i in seq_along(mapbiomas_classes_to_keep)) {
    class <- mapbiomas_classes_to_keep[i]
    window_recoded[window == class] <- i
}

## Rare classes are recoded as NA
message("Recoding the following rare classes as NA:")
print(rare_mapbiomas_classes)
window_recoded[window %in% rare_mapbiomas_classes] <- NA

table(window)
table(window_recoded)
mean(is.na(window_recoded))

message("Fraction of pixels with at least one missing value in the recoded data:")
mean(rowMeans(is.na(window_recoded)) > 0)
message("Fraction of pixels missing in >50% of years in the recoded data:")
mean(rowMeans(is.na(window_recoded)) > .5)
message("Fraction of pixels missing in 100% of years in the recoded data:")
mean(rowMeans(is.na(window_recoded)) == 1.0)

full_panel <- apply(window_recoded, 1, function(y) list(y=as.vector(y), time=seq_along(y)))

## Using prob=valid_pixel_index excludes pixels that have 100% missing observations in the original data
panel <- sample(full_panel, size=length(full_panel) * opt$subsample, replace=FALSE, prob=valid_pixel_index)

## We no longer need the full window at this point, rm it to save memory
rm(window)
rm(full_panel)
gc()

## These aren't actually used in optimization,
## they're just used to create other parameters of the same shape/dimension/time horizon
n_states <- length(mapbiomas_classes_to_keep)
n_time_periods <- ncol(window_recoded)
dummy_pr_transition <- 0.2 * matrix(1/n_states, nrow=n_states, ncol=n_states) + 0.8 * diag(n_states)
dummy_pr_y <- 0.2 * matrix(1/n_states, n_states, n_states) + 0.8 * diag(n_states)
dummy_params <- list(mu=rep(1/n_states, n_states),
                     P_list=rep(list(dummy_pr_transition), n_time_periods - 1),
                     pr_y=dummy_pr_y,
                     n_components=n_states)

estimates <- get_em_and_min_dist_estimates_random_initialization(params=dummy_params,
                                                                 panel=panel,
                                                                 n_random_starts_em=opt$n_random_starts_em,
                                                                 n_random_starts_md=opt$n_random_starts_md,
                                                                 diag_min=0.8,
                                                                 diag_max=0.95,
                                                                 skip_ml_if_md_is_diag_dominant=opt$skip_ml_if_md_is_diag_dominant,
                                                                 use_md_as_initial_values_for_em=opt$use_md_as_initial_values_for_em)

estimates$P_hat_frequency <- lapply(estimates$M_Y_joint_hat, get_transition_probs_from_M_S_joint)

estimates$mapbiomas_classes_to_keep <- mapbiomas_classes_to_keep
estimates$rare_mapbiomas_classes <- rare_mapbiomas_classes
estimates$class_frequencies <- class_frequencies
estimates$class_frequencies_before_combining <- class_frequencies_before_combining

estimates$options <- opt
estimates$window_bbox <- as.data.frame(bbox(window_extent))

estimates$fraction_missing_in_all_years <- fraction_missing_in_all_years
estimates$count_missing_in_all_years <- count_missing_in_all_years

filename <- sprintf("./atlantic_forest_output/estimates_window_%s_%s_width_%s_class_frequency_cutoff_%s_subsample_%s_combined_classes%s%s%s%s.rds",
                    opt$row, opt$col, opt$width_in_pixels, opt$class_frequency_cutoff, opt$subsample,
                    ifelse(opt$grassland_as_forest, "_grassland_as_forest", ""),
                    ifelse(opt$combine_other_non_forest, "_combine_other_non_forest", ""),
                    ifelse(opt$skip_ml_if_md_is_diag_dominant, "_skip_ml_if_md_is_diag_dominant", ""),
                    ifelse(opt$use_md_as_initial_values_for_em, "_use_md_as_initial_values_for_em", ""))
message("Saving ", filename)
saveRDS(estimates, file=filename)

for(class_index in seq_along(estimates$mapbiomas_classes_to_keep)) {
    class <- estimates$mapbiomas_classes_to_keep[class_index]

    ## Diagonals of the transition matrix (for example, Pr[ forest at t+1 | forest at t ])
    P_hat_frequency <- sapply(estimates$P_hat_frequency, function(P) P[class_index, class_index])
    P_hat_md <- sapply(estimates$min_dist_params_hat_best_objfn$P_list, function(P) P[class_index, class_index])

    if("em_params_hat_best_likelihood" %in% names(estimates)) {
        P_hat_ml <- sapply(estimates$em_params_hat_best_likelihood$P_list, function(P) P[class_index, class_index])
    } else {
        P_hat_ml <- rep(NA, length(P_hat_md))
    }

    df <- data.table(time_index=seq_along(P_hat_frequency), P_hat_frequency=P_hat_frequency, P_hat_md, P_hat_ml)
    df_melted <- melt(df, id.vars="time_index")

    title <- sprintf("Probability of Remaining in Mapbiomas Class %s", class)  # TODO Window info in title?
    p <- (ggplot(df_melted, aes(x=time_index, y=value, group=variable, color=variable)) +
          geom_point() +
          geom_line() +
          ggtitle(title) +
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_color_discrete("algorithm") +
          ylab("probability") +
          theme_bw())
    filename <- sprintf("transition_matrix_diagonals_window_%s_%s_width_%s_class_%s_with_combined_classes_%s.png",
                        opt$row, opt$col, opt$width_in_pixels, class, ifelse(opt$grassland_as_forest, "grassland_as_forest", ""))
    ggsave(p, filename=filename, width=6, height=4, units="in")
}
