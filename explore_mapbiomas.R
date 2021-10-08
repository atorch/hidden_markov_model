library(Rsolnp)
library(data.table)
library(optparse)
library(raster)

source("src/hmm_functions.R")

mapbiomas <- stack("HMM_MapBiomas_v2/mapbiomas.vrt")

## TODO Run this script in parallel with different input options
opt_list <- list(make_option("--row", default=90400, type="integer"),
                 make_option("--col", default=23500, type="integer"),
                 make_option("--width_in_pixels", default=200, type="integer"),
                 make_option("--class_frequency_cutoff", default=0.005, type="double"))
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

## Note that the set of unique observations varies by year!
sort(unique(window[, 3]))
sort(unique(window[, 4]))

## Some observations are extremely rare (11 and 12 and 33 for example)
table(window[, 4])
table(window[, 4]) / nrow(window)

table(window)
round(table(window) / (nrow(window) * ncol(window)), 4)

## Careful, there can be missing values!
mean(is.na(c(window, recursive=TRUE)))
sum(is.na(c(window, recursive=TRUE)))

window_recoded <- window

## TODO We could collapse multiple mapbiomas classes together (instead of replacing rare classes with NAs)

for(i in seq_along(mapbiomas_classes_to_keep)) {
    class <- mapbiomas_classes_to_keep[i]
    window_recoded[window == class] <- i
}

## Rare classes are recoded as NA
window_recoded[window %in% rare_mapbiomas_classes] <- NA

table(window)
table(window_recoded)
mean(is.na(window_recoded))

panel <- apply(window_recoded, 1, function(y) list(y=as.vector(y), time=seq_along(y)))

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

for(idx in seq_along(panel)) {
    panel[[idx]]$point_id <- idx
    panel[[idx]]$time <- seq_along(panel[[idx]]$y)
}

dtable <- rbindlist(Map(data.frame, panel))
setkey(dtable, point_id)

stopifnot(all(c("point_id", "time", "y") %in% names(dtable)))

dtable[, y_one_period_ahead := c(tail(y, .N-1), NA), by="point_id"]
dtable[, y_two_periods_ahead := c(tail(y, .N-2), NA, NA), by="point_id"]

## Joint distribution of (Y_{t+1}, Y_{t})
M_Y_joint_hat_list <- lapply(seq_len(max(dtable$time) - 1), function(fixed_t) {
    with(subset(dtable, time == fixed_t), prop.table(table(y_one_period_ahead, y)))
})

## Compute inverses once and pass them to get_min_distance_estimates / solnp
M_Y_joint_hat_inverse_list <- lapply(M_Y_joint_hat_list, solve)

## Joint distribution of (Y_{t+2}, Y_{t+1}, Y_{t})
M_fixed_y_Y_joint_hat_list <- lapply(seq_len(dummy_params$n_components), function(fixed_y) {
    lapply(seq_len(max(dtable$time) - 2), function(fixed_t) {
        ## Note: we need to pass factors to table() so that it includes
        ## rows and columns of zeros in cases where a certain class (factor level) isn't observed
        levels <- seq_len(dummy_params$n_components)
        return(with(subset(dtable, time == fixed_t & y_two_periods_ahead == fixed_y),
                    table(factor(y_one_period_ahead, levels=levels),
                          factor(y, levels=levels))) / sum(dtable$time == fixed_t &
                                                           !is.na(dtable$y_two_periods_ahead) &
                                                           !is.na(dtable$y_one_period_ahead) &
                                                           !is.na(dtable$y)))
    })
})

md_estimates <- get_min_distance_estimates(dummy_params, M_Y_joint_hat_list, M_Y_joint_hat_inverse_list, M_fixed_y_Y_joint_hat_list, dtable)

## TODO Do we really need 10 initializations for EM?  Do they end up in different places?  Make it a fn argument (and optparse)
estimates <- get_hmm_and_minimum_distance_estimates_random_initialization(params=dummy_params, panel=panel, n_random_starts=5)

estimates$em_params_hat_best_likelihood

estimates$mapbiomas_classes_to_keep <- mapbiomas_classes_to_keep
estimates$rare_mapbiomas_classes <- rare_mapbiomas_classes

estimates$P_hat_frequency <- lapply(estimates$M_Y_joint_hat, get_transition_probs_from_M_S_joint)

filename <- sprintf("estimates_window_%s_%s_width_%s_class_frequency_cutoff_%s.rds",
                    opt$row, opt$col, opt$width_in_pixels, opt$class_frequency_cutoff)
saveRDS(estimates, file=filename)
