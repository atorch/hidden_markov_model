library(Rsolnp)
library(data.table)
library(ggplot2)
library(optparse)
library(raster)

source("src/hmm_functions.R")

mapbiomas <- stack("HMM_MapBiomas_v2/mapbiomas.vrt")

## TODO Run this script in parallel with different input options
opt_list <- list(make_option("--row", default=90000, type="integer"),
                 make_option("--col", default=23000, type="integer"),
                 make_option("--width_in_pixels", default=500, type="integer"),
                 make_option("--subsample", default=0.1, type="double"),
                 make_option("--class_frequency_cutoff", default=0.005, type="double"),
                 make_option("--n_random_starts", default=5, type="integer"))
opt <- parse_args(OptionParser(option_list=opt_list))
message("command line options: ", paste(sprintf("%s=%s", names(opt), opt), collapse=", "))

nlayers(mapbiomas)

window <- getValuesBlock(mapbiomas,
                         row=opt$row,
                         col=opt$col,
                         nrows=opt$width_in_pixels,
                         ncols=opt$width_in_pixels)
dim(window)

window_extent <- extent(mapbiomas, opt$row, opt$row + opt$width_in_pixels, opt$col, opt$col + opt$width_in_pixels)
window_raster<- raster(window_extent, crs=crs(mapbiomas), nrows=opt$width_in_pixels, ncols=opt$width_in_pixels)
values(window_raster) <- window[, 1]

filename <- sprintf("raster_window_%s_%s_width_%s_band_1.tif", opt$row, opt$col, opt$width_in_pixels)
writeRaster(window_raster, filename, overwrite=TRUE)

## Unique land use classes in this window are 3, 9, 11, 12, 21, 22, 29, 33
## See https://mapbiomas-br-site.s3.amazonaws.com/downloads/Colecction%206/Cod_Class_legenda_Col6_MapBiomas_BR.pdf
unique_mapbiomas_classes <- sort(unique(c(window, recursive=TRUE)))

rare_mapbiomas_classes <- vector("numeric")
for(class in unique_mapbiomas_classes) {
    if(mean(window == class, na.rm=TRUE) < opt$class_frequency_cutoff) {
        rare_mapbiomas_classes <- c(rare_mapbiomas_classes, class)
    }
}

## TODO We are going to recode rare classes as NA
## This is effectively assuming that all observations of rare classes must be misclassifications
## Keeps classes 3 and 21 (forest and mosaic of pasture + agriculture)
mapbiomas_classes_to_keep <- unique_mapbiomas_classes[!unique_mapbiomas_classes %in% rare_mapbiomas_classes]

## TODO Drop points that are NA in every year
## TODO Count the number of points that are NA in every year
## TODO Sanity check, run Viterbi, save rasters, visualize

## Note that the set of unique observations varies by year!
sort(unique(window[, 3]))
sort(unique(window[, 4]))

## Some observations are extremely rare (11 and 12 and 33 for example)
table(window[, 4])
table(window[, 4]) / nrow(window)

table(window)  # TODO Logic for skipping windows that aren't mainly forest + crops (skip if too much grassland, urban, etc)
round(table(window) / (nrow(window) * ncol(window)), 4)

## Careful, there can be missing values!
mean(is.na(c(window, recursive=TRUE)))
sum(is.na(c(window, recursive=TRUE)))

window_recoded <- window

for(i in seq_along(mapbiomas_classes_to_keep)) {
    class <- mapbiomas_classes_to_keep[i]
    window_recoded[window == class] <- i
}

## Rare classes are recoded as NA
window_recoded[window %in% rare_mapbiomas_classes] <- NA

table(window)
table(window_recoded)
mean(is.na(window_recoded))

full_panel <- apply(window_recoded, 1, function(y) list(y=as.vector(y), time=seq_along(y)))

panel <- sample(full_panel, size=length(full_panel) * opt$subsample, replace=FALSE)

## These aren't actually used in optimization,
## they're just used to create other parameters of the same shape/dimension/time horizon
n_states <- length(mapbiomas_classes_to_keep)
n_time_periods <- ncol(window)
dummy_pr_transition <- 0.2 * matrix(1/n_states, nrow=n_states, ncol=n_states) + 0.8 * diag(n_states)
dummy_pr_y <- 0.2*matrix(1/n_states, n_states, n_states) + 0.8*diag(n_states)
dummy_params <- list(mu=rep(1/n_states, n_states),
                     P_list=rep(list(dummy_pr_transition), n_time_periods - 1),
                     pr_y=dummy_pr_y,
                     n_components=n_states)

estimates <- get_hmm_and_minimum_distance_estimates_random_initialization(params=dummy_params,
                                                                          panel=panel,
                                                                          n_random_starts=opt$n_random_starts)

estimates$em_params_hat_best_likelihood

estimates$mapbiomas_classes_to_keep <- mapbiomas_classes_to_keep
estimates$rare_mapbiomas_classes <- rare_mapbiomas_classes

estimates$P_hat_frequency <- lapply(estimates$M_Y_joint_hat, get_transition_probs_from_M_S_joint)  # TODO Also save mapbiomas class counts (tables) in estimates

estimates$options <- opt

filename <- sprintf("estimates_window_%s_%s_width_%s_class_frequency_cutoff_%s_subsample_%s.rds",
                    opt$row, opt$col, opt$width_in_pixels, opt$class_frequency_cutoff, opt$subsample)
saveRDS(estimates, file=filename)

for(class_index in seq_along(estimates$mapbiomas_classes_to_keep)) {
    class <- estimates$mapbiomas_classes_to_keep[class_index]

    ## Diagonals of the transition matrix (for example, Pr[ forest at t+1 | forest at t ])
    P_hat_frequency <- sapply(estimates$P_hat_frequency, function(P) P[class_index, class_index])
    P_hat_md <- sapply(estimates$min_dist_params_hat_best_objfn$P_list, function(P) P[class_index, class_index])
    P_hat_ml <- sapply(estimates$em_params_hat_best_likelihood$P_list, function(P) P[class_index, class_index])

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
    filename <- sprintf("transition_matrix_diagonals_window_%s_%s_width_%s_class_%s.png", opt$row, opt$col, opt$width_in_pixels, class)
    ggsave(p, filename=filename, width=6, height=4, units="in")
}
