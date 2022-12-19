#!/usr/bin/env Rscript

rm(list =ls())
library(parallel)
library(raster)
library(data.table)
library(stringr)
library(optparse)

source("hmm_functions.R")

opt_list <- list(make_option("--mapbiomas_raster_path", default="./HMM_MapBiomas_v2/mapbiomas.vrt"),
                 make_option("--row", default=54001, type="integer"),
                 make_option("--col", default=5001, type="integer"),
                 make_option("--width_in_pixels", default=1000, type="integer"),
                 make_option("--raster_year", default=2017, type ="integer"),
                 make_option("--subsample", default=0.01, type="double"))

opt <- parse_args(OptionParser(option_list=opt_list))
message("command line options: ", paste(sprintf("%s=%s", names(opt), opt), collapse=", "))

yearsVec <- 1985:2020
yearsVecExpand <- 1985:2040
row <- opt$row
col <- opt$col
width_in_pixels <- opt$width_in_pixels
subsample_for_viterbi <- opt$subsample
rasterYear <- opt$raster_year

## This step uses .rds files saved by run_estimation_single_mapbiomas_window.R
## These files are not included in the repo
filename <- sprintf("estimates_window_%s_%s_width_%s_class_frequency_cutoff_0.005_subsample_0.01_combined_classes_grassland_as_forest_combine_other_non_forest_use_md_as_initial_values_for_em.rds", row, col, width_in_pixels)
fullFilePath <- file.path("atlantic_forest_output", filename)
if (!file.exists(fullFilePath)){
    message("No Estimates for this Window -- exiting")
    q("no")
}
estimates <- readRDS(fullFilePath)

## Load mapbiomas raster
mapbiomas <- stack(opt$mapbiomas_raster_path)

## Carbon Stock File (file is Carbon Stock in rasterYear)
## This file is not included in the repo
carbonFile <- paste0('./carbon_stock_data/carbonStockRaster', rasterYear, '.tif')
carbonRaster <- terra::rast(carbonFile)

## Stop analysis if pr_y is not diagonally dominant
if(any(diag(estimates$em_params_hat_best_likelihood$pr_y) < .5)) {
    message("The diagonals of estimates$em_params_hat_best_likelihood$pr_y are ",
            paste(diag(estimates$em_params_hat_best_likelihood$pr_y), collapse=", "))
     q("no")
}

## Get values for analysis from mapbiomas raster
window <- getValuesBlock(mapbiomas,
                         row=row,
                         col=col,
                         nrows=width_in_pixels,
                         ncols=width_in_pixels)

cellsToPull <- cellFromRowColCombine(mapbiomas, seq(row,row+width_in_pixels-1),seq(col,col+width_in_pixels-1))

opt <- estimates$opt
mapbiomas_classes_to_keep <- estimates$mapbiomas_classes_to_keep

## We need to combine classes (in the exact same way we did before estimation)
## and then generate the recoded window before running Viterbi
## TODO A bunch of this is copy pasted
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
pixelsToUse <- sample.int(length(full_panel), size=length(full_panel) * subsample_for_viterbi, replace=FALSE, prob=valid_pixel_index)
panel <- full_panel[pixelsToUse]

##Viterbi
viterbi <- apply_viterbi_path_in_parallel(panel, params_hat=estimates$em_params_hat_best_likelihood, max_cores=1)

## Recode as as mapbiomas classes, use estimates$mapbiomas_classes_to_keep
forest_class <- 3
agriculture_and_pasture_class <- 21
water_rocks_class <- 33

## Simulate path forward for another 20 years
transMatrix <- estimates$em_params_hat_best_likelihood$P_list[[35]]
simMarkovForwardFunc <- function(currentState, x) {
    sample.int(n=estimates$em_params_hat_best_likelihood$n_components, size=1, prob=transMatrix[currentState,])
}

nYearSim <- yearsVecExpand[length(yearsVecExpand)]-yearsVec[length(yearsVec)]
futureSimBase <- lapply(viterbi,
                        function(z) Reduce(simMarkovForwardFunc, x = 1:nYearSim, init = z[length(z)],accumulate = TRUE))

##counterfactual where forest is an absorbing state
forestIndex <- which(estimates$mapbiomas_classes_to_keep == forest_class)
transMatrix[1,1] <- 1
transMatrix[1,2:length(estimates$mapbiomas_classes_to_keep)] <-0 

##Process Viterbi Output
futureSimNoDefor <- lapply(viterbi,
                        function(z) Reduce(simMarkovForwardFunc, x = 1:nYearSim, init = z[length(z)],accumulate = TRUE))


## pixelTmpIndex is the row number of the pixel in the panel list
yearVecTmp <- yearsVecExpand[(length(yearsVec)):length(yearsVecExpand)]
processedDatWideList <- list(
    viterbiDT = as.data.table(viterbi)[,year := yearsVec],
    futureSimBaseDT = as.data.table(futureSimBase)[,year := yearVecTmp],
    futureSimNoDeforDT = as.data.table(futureSimNoDefor)[,year := yearVecTmp],
    obsDT = as.data.table(lapply(panel,function(dat) dat$y))[,year := yearsVec])

processedDatDTList <- lapply(processedDatWideList,
                         function(x)
                             melt(x,
                                  id.var = 'year',
                                  variable.name = 'pixelTmpIndex',
                                  value.name = 'recodedLandUse')
                         )

setnames(processedDatDTList$futureSimBaseDT,'recodedLandUse','recodedLandUse.viterbi_base')
setnames(processedDatDTList$futureSimNoDeforDT,'recodedLandUse','recodedLandUse.viterbi_nodefor')
setnames(processedDatDTList$viterbiDT,'recodedLandUse','recodedLandUse.viterbi')
setnames(processedDatDTList$obsDT,'recodedLandUse','recodedLandUse.y')

landUseOverTime <-
    Reduce(function(dat,x) merge(dat,x, by = c('year','pixelTmpIndex'),all=TRUE),processedDatDTList)

landUseOverTime[year%in% yearsVec,':='(recodedLandUse.viterbi_base = recodedLandUse.viterbi,
                                       recodedLandUse.viterbi_nodefor = recodedLandUse.viterbi)]

landUseOverTime[,pixelTmpIndex  := as.integer(str_replace(pixelTmpIndex,'V',''))]

##Map from pixelTmpIndex to pixelOrigIndex (from index in calculations mapbiomas raster)
pixelMapping <- data.table(pixelTmpIndex = 1:length(pixelsToUse),
                           pixelWindowIndex = pixelsToUse,
                           pixelMapBiomasIndex = cellsToPull[pixelsToUse])
landUseOverTime <- merge(landUseOverTime, pixelMapping, by = 'pixelTmpIndex',all.x=TRUE)



landUseLevels <- c(forest = which(estimates$mapbiomas_classes_to_keep == forest_class),
                   ag_past = which(estimates$mapbiomas_classes_to_keep == agriculture_and_pasture_class),
                   water_rocks = which(estimates$mapbiomas_classes_to_keep == agriculture_and_pasture_class))
landUseOverTime[,recodedLandUse.viterbiF := factor(recodedLandUse.viterbi,levels = landUseLevels,labels = names(landUseLevels))]
landUseOverTime[,recodedLandUse.yF := factor(recodedLandUse.y,levels = landUseLevels,labels = names(landUseLevels))]
landUseOverTime[,recodedLandUse.viterbi_baseF := factor(recodedLandUse.viterbi_base,levels = landUseLevels,labels = names(landUseLevels))]
landUseOverTime[,recodedLandUse.viterbi_nodeforF := factor(recodedLandUse.viterbi_nodefor,levels = landUseLevels,labels = names(landUseLevels))]


##Recover estimates of carbon stock for different age forests (using year rasterYear map of carbon stock)
landUseOverTime[,forestDumV := recodedLandUse.viterbiF == 'forest']
landUseOverTime[,forestDumY := recodedLandUse.yF == 'forest']
landUseOverTime[,forestDumVB := recodedLandUse.viterbi_baseF == 'forest']
landUseOverTime[,forestDumVND := recodedLandUse.viterbi_nodeforF == 'forest']

laggedValsV <- landUseOverTime[,data.table::shift(.SD,0:(max(year)-min(year)),give.names=TRUE),.SDcols = 'forestDumV',by = 'pixelMapBiomasIndex']
laggedValsY <- landUseOverTime[,data.table::shift(.SD,0:(max(year)-min(year)),give.names=TRUE),.SDcols = 'forestDumY',by = 'pixelMapBiomasIndex']
laggedValsVB <- landUseOverTime[,data.table::shift(.SD,0:(max(year)-min(year)),give.names=TRUE),.SDcols = 'forestDumVB',by = 'pixelMapBiomasIndex']
laggedValsVND <- landUseOverTime[,data.table::shift(.SD,0:(max(year)-min(year)),give.names=TRUE),.SDcols = 'forestDumVND',by = 'pixelMapBiomasIndex']

setnames(laggedValsY,'pixelMapBiomasIndex','pixelMapBiomasIndexCheck')


##Get age of forest
forestAgeDat <- cbind(laggedValsV,laggedValsY,laggedValsVB, laggedValsVND, landUseOverTime[,list(year)])


##The Position function gives the first lag in which something other than forest is in the data (so that is the age of the forest)
forestAgeDat[forestDumV_lag_0 == TRUE, forest_ageV := apply(.SD, MARGIN = 1, FUN = Position,f=function(x) ifelse(is.na(x),FALSE,!x), nomatch = Inf),
             .SDcols = grep('forestDumV_lag_[1-9]{1}[0-9]{0,1}',names(forestAgeDat))]
forestAgeDat[forestDumY_lag_0 == TRUE, forest_ageY := apply(.SD, MARGIN = 1, FUN = Position,f=function(x) ifelse(is.na(x),FALSE,!x), nomatch = Inf),
               .SDcols = grep('forestDumY_lag_[1-9]{1}[0-9]{0,1}',names(forestAgeDat))]
forestAgeDat[forestDumVB_lag_0 == TRUE, forest_ageVB := apply(.SD, MARGIN = 1, FUN = Position,f=function(x) ifelse(is.na(x),FALSE,!x), nomatch = Inf),
             .SDcols = grep('forestDumVB_lag_[1-9]{1}[0-9]{0,1}',names(forestAgeDat))]
forestAgeDat[forestDumVND_lag_0 == TRUE, forest_ageVND := apply(.SD, MARGIN = 1, FUN = Position,f=function(x) ifelse(is.na(x),FALSE,!x), nomatch = Inf),
               .SDcols = grep('forestDumVND_lag_[1-9]{1}[0-9]{0,1}',names(forestAgeDat))]


set(forestAgeDat, ,grep('forestDumV_lag_[1-9]{1}[0-9]{0,1}',names(forestAgeDat)),NULL)
set(forestAgeDat, ,grep('forestDumY_lag_[1-9]{1}[0-9]{0,1}',names(forestAgeDat)),NULL)
set(forestAgeDat, ,grep('forestDumVND_lag_[1-9]{1}[0-9]{0,1}',names(forestAgeDat)),NULL)
set(forestAgeDat, ,grep('forestDumVB_lag_[1-9]{1}[0-9]{0,1}',names(forestAgeDat)),NULL)



##Get carbon stock for landuses
##first forest
carbonStockDatForest <- forestAgeDat[year == rasterYear & forestDumV_lag_0 == TRUE]
carbonStockDatForest[,c('xCoord', 'yCoord') := data.frame(xyFromCell(mapbiomas, pixelMapBiomasIndex))]
carbonStockDatForest[,carbonVal := terra::extract(carbonRaster, cbind(x=xCoord,y=yCoord))]

##Get CarbonStock for non-forest
carbonStockDatNonForest <- landUseOverTime[(forestDumV == FALSE |is.na(forestDumV))& year == rasterYear,
                                           list(pixelMapBiomasIndex,recodedLandUse.viterbiF)]
carbonStockDatNonForest[,c('xCoord', 'yCoord') := data.frame(xyFromCell(mapbiomas, pixelMapBiomasIndex))]
carbonStockDatNonForest[,carbonVal := terra::extract(carbonRaster, cbind(x=xCoord,y=yCoord))]


##Get averages of carbon stock for forest and non-forest
carbonStockDatNonForest[,list(avg = mean(carbonVal,na.rm=TRUE),
                              sd = sd(carbonVal,na.rm = TRUE)),
                        by = recodedLandUse.viterbiF]

carbonStockDatForest[,list(avg = mean(carbonVal,na.rm=TRUE),
                           sd = sd(carbonVal,na.rm=TRUE))]

viterbiResults <- list(csNonForest = carbonStockDatNonForest,csForest = carbonStockDatForest,landuse = forestAgeDat)

saveRDS(viterbiResults, file.path("carbon_stock_results", sprintf("landUseAndCarbon_%s_%s_width_%s_%s.rds", row, col, width_in_pixels, rasterYear)))
