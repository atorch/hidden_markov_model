library(parallel)
library(raster)

source("src/hmm_functions.R")

row <- 76001
col <- 28001
width_in_pixels <- 1000

filename <- sprintf("estimates_window_%s_%s_width_%s_class_frequency_cutoff_0.005_subsample_0.01_combined_classes_grassland_as_forest_combine_other_non_forest_use_md_as_initial_values_for_em.rds", row, col, width_in_pixels)

estimates <- readRDS(filename)

subsample_for_viterbi <- 0.01

mapBioMassFile <- "./HMM_MapBiomas_v2/mapbiomas.vrt"
mapbiomas <- stack(mapBioMassFile)

window <- getValuesBlock(mapbiomas,
                         row=row,
                         col=col,
                         nrows=width_in_pixels,
                         ncols=width_in_pixels)

dim(window)

opt <- estimates$opt
mapbiomas_classes_to_keep <- estimates$mapbiomas_classes_to_keep

## We need to combine classes (in the exact same way we did before estimation)
## and then generate the recoded window before running Viterbi
## A bunch of this is copy pasted from explore_mapbiomas.R
## TODO Put it in a function that we can import/source and reuse

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

## Note: the code expects observations to be in the set {1, 2, 3, ..., |Y|},
## So we need to recode sets of classes like {3, 21} to {1, 2} for example
window_recoded <- window

for(i in seq_along(mapbiomas_classes_to_keep)) {
    class <- mapbiomas_classes_to_keep[i]
    window_recoded[window == class] <- i
}

## Rare classes are recoded as NA
window_recoded[window %in% rare_mapbiomas_classes] <- NA

full_panel <- apply(window_recoded, 1, function(y) list(y=as.vector(y), time=seq_along(y)))

## Using prob=valid_pixel_index excludes pixels that have 100% missing observations in the original data
panel <- sample(full_panel, size=length(full_panel) * subsample_for_viterbi, replace=FALSE, prob=valid_pixel_index)

viterbi <- apply_viterbi_path_in_parallel(panel, params_hat=estimates$em_params_hat_best_likelihood, max_cores=30)

## This is the most likely sequence of hidden states for the first pixel
## Compare to panel[[1]]$y
viterbi[[1]]

## To see these as mapbiomas classes, use estimates$mapbiomas_classes_to_keep
## For example:
forest_class <- 3
agriculture_and_pasture_class <- 21
forest_index <- which(estimates$mapbiomas_classes_to_keep == forest_class)

## When was the first pixel forest?
viterbi[[1]] == forest_index
