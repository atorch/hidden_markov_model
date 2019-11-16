library(data.table)
library(digest)
library(gbm)
library(ggplot2)
library(optparse)
library(plyr); library(dplyr)  # For bind_rows
library(randomForest)
library(Rsolnp)
library(sp)

set.seed(789)

# Binary for {crops, pasture}, ternary for {soy, other crops, pasture}
opt_list <- list(make_option("--landuse_set", default="binary"),
                 make_option("--panel_length", default="short"),
                 make_option("--classifier_training_fraction", default=0.1, type="double"),
                 make_option("--classifier_pasture_fraction", default=0.6, type="double"))
opt <- parse_args(OptionParser(option_list=opt_list))
message("command line options: ", paste(sprintf("%s=%s", names(opt), opt), collapse=", "))

stopifnot(opt$landuse_set %in% c("binary", "ternary"))
stopifnot(opt$panel_length %in% c("long", "short"))
stopifnot(0 < opt$classifier_training_fraction && opt$classifier_training_fraction < 1)

source("embrapa_validation_landuse_mapping.R")
source("hmm_functions.R")
source("ggplot_utils.R")
set_ggplot_theme()

points <- readRDS("~/Dropbox/amazon_hmm_shared/embrapa_validation/embrapa_validation_points_with_covariates.rds")
message("done loading ", nrow(points), " point-years, ", length(unique(points$point_id)), " unique spatial points")

n_points_before_removing <- length(unique(points$point_id))

points_with_validation_landuse <- readRDS("~/Dropbox/amazon_hmm_shared/embrapa_validation/embrapa_validation_points_with_validation_landuse.rds")
points <- merge(points, points_with_validation_landuse, all.x=TRUE, by=c("point_id", "year"))  # Adds validation_landuse column
points <- arrange(points, point_id, year)

state_abbr <- unique(as.character(points$state))  # E.g. MT for Mato Grosso
stopifnot(length(state_abbr) == 1)
stopifnot(nchar(state_abbr) == 2)
message("all points are in ", state_abbr)

validation_years <- sort(unique(points$year[!is.na(points$validation_landuse)]))
message("years with validation land use data: ", paste(validation_years, collapse=", "))

if(opt$panel_length == "long") {
    ## More precise HMM parameter estimates than short panel
    points <- subset(points, year >= min(validation_years))
    
} else {
    points <- subset(points, year %in% validation_years)
}
message("keeping the following years: ", paste(unique(points$year), collapse=", "))
message("keeping ", nrow(points), " point-years, ", length(unique(points$point_id)), " unique spatial points")

point_id_ever_mata <- unique(points$point_id[!is.na(points$validation_landuse) & tolower(points$validation_landuse) == "mata"])
point_id_ever_reflorestamento <- unique(points$point_id[!is.na(points$validation_landuse) & tolower(points$validation_landuse) == "reflorestamento"])
points <- subset(points, !point_id %in% c(point_id_ever_mata, point_id_ever_reflorestamento))  # Down to 403 points
message("keeping ", nrow(points), " point-years, ", length(unique(points$point_id)), " unique spatial points after removing points whose validation landuse was ever mata or reflorestamento")

message("running validation exercise with ", opt$landuse_set, " land use set S")
if(opt$landuse_set == "binary") {
    map_landuse_to_S <- map_landuse_to_S_binary
    points$validation_landuse_coarse <- plyr::revalue(points$validation_landuse, replace=map_landuse_to_S)
    stopifnot(all(is.na(points$validation_landuse_coarse) | points$validation_landuse_coarse %in% c("crops", "pasture")))
} else {
    map_landuse_to_S <- map_landuse_to_S_ternary
    points$validation_landuse_coarse <- plyr::revalue(points$validation_landuse, replace=map_landuse_to_S)
    stopifnot(all(is.na(points$validation_landuse_coarse) | points$validation_landuse_coarse %in% c("soy", "non-soy crops", "pasture")))
}

table(points$validation_landuse_coarse)  # Careful, pasture is rare
points$validation_landuse_coarse <- factor(points$validation_landuse_coarse)  # Factor for randomForest classification

## Note: pasture would be very rare in the training set if we took a simple random sample,
## so we stratify in order to over-sample pasture
n_point_id_train <- floor(opt$classifier_training_fraction * length(unique(points$point_id)))
subset_pasture <- subset(points, validation_landuse_coarse == "pasture")
point_id_train_pasture <- sample(unique(subset_pasture$point_id), size=n_point_id_train * opt$classifier_pasture_fraction)
point_id_train <- c(point_id_train_pasture, sample(subset(points, !point_id %in% point_id_train_pasture)$point_id, size=n_point_id_train - length(point_id_train_pasture)))

## MODIS variables
evi_vars <- sprintf("evi_%s", seq(1, 23))
nir_vars <- sprintf("nir_%s", seq(1, 23))
mir_vars <- sprintf("mir_%s", seq(1, 23))
red_vars <- sprintf("red_%s", seq(1, 23))
blue_vars <- sprintf("blue_%s", seq(1, 23))
stopifnot(all(evi_vars %in% names(points)))
stopifnot(all(nir_vars %in% names(points)))
stopifnot(all(mir_vars %in% names(points)))
stopifnot(all(red_vars %in% names(points)))
stopifnot(all(blue_vars %in% names(points)))

get_initial_hmm_params <- function(landuse_set) {
    if(landuse_set == "binary") {
        initial_hmm_params <- list(list(P=rbind(c(0.90, 0.10),
                                                c(0.10, 0.90)),
                                        mu=c(0.6, 0.4),
                                        pr_y=rbind(c(0.80, 0.20),
                                                   c(0.20, 0.80)),
                                        n_components=2),
                                   list(P=rbind(c(0.90, 0.10),
                                                c(0.30, 0.70)),
                                        mu=c(0.6, 0.4),
                                        pr_y=rbind(c(0.90, 0.10),
                                                   c(0.30, 0.70)),
                                        n_components=2),
                                   list(P=rbind(c(0.90, 0.10),
                                                c(0.40, 0.60)),
                                        mu=c(0.6, 0.4),
                                        pr_y=rbind(c(0.90, 0.10),
                                                   c(0.30, 0.70)),
                                        n_components=2))
    } else {
        initial_hmm_params <- list(list(P=rbind(c(0.80, 0.10, 0.10),
                                                c(0.10, 0.80, 0.10),
                                                c(0.10, 0.10, 0.80)),
                                        mu=c(1/3, 1/3, 1/3),
                                        pr_y=rbind(c(0.80, 0.10, 0.10),
                                                   c(0.10, 0.80, 0.10),
                                                   c(0.10, 0.10, 0.80)),
                                        n_components=3),
                                   list(P=rbind(c(.60, 0.20, 0.20),
                                                c(0.20, 0.60, 0.20),
                                                c(0.10, 0.10, 0.80)),
                                        mu=c(0.2, 0.3, 0.5),
                                        pr_y=rbind(c(0.60, 0.20, 0.20),
                                                   c(0.20, 0.60, 0.20),
                                                   c(0.10, 0.10, 0.80)),
                                        n_components=3))
    }
    return(initial_hmm_params)
}

points_train <- subset(points, point_id %in% point_id_train)
points_test <- subset(points, !point_id %in% point_id_train)

table(points_train$validation_landuse_coarse)  # Careful, pasture is rare -- should make sure that training set is not 100% crops, same for test
table(points_test$validation_landuse_coarse)  # Similar comment about non-soy crops in ternary case

predictors <- c(evi_vars, nir_vars, mir_vars, red_vars, blue_vars)
model_formula <- reformulate(termlabels=predictors, response="validation_landuse_coarse")

## Fit GBM (tried random forest but confusion matrix was generally not diagonally dominant, which violates the HMM assumptions)
gbm_training_set <- subset(points_train, !is.na(validation_landuse_coarse))
message("training ML model (GBM classifier) on ", length(unique(gbm_training_set$point_id)), " points")
model_gbm <- gbm(formula=model_formula, distribution="multinomial",  # Does this have to be bernoulli for two-class case?
                 data=gbm_training_set, n.trees=5000, interaction.depth=4, shrinkage=0.001, cv.folds=3)

ntrees_star <- gbm.perf(model_gbm, method="cv")  # Plot shows training and CV deviance as function of number of trees -- selects around 5000-6000 trees
message("GBM uses ", ntrees_star, " trees (tuning parameter selected by cross-validation)")

gbm_predictions <- predict(model_gbm, type="response", newdata=points_test, ntrees=ntrees_star)  # Outputs class probabilities
points_test$landuse_predicted_gbm <- colnames(gbm_predictions)[as.integer(apply(gbm_predictions, 1, which.max))]
mean(is.na(points_test$landuse_predicted_gbm))  # Zero (i.e. no missing predictions), unlike random forest

table(points_test$landuse_predicted_gbm)

## Careful, not diagonally dominant in ternary case (see non-soy crops entry), expect HMM to perform poorly
gbm_test_confusion <- table(points_test$validation_landuse_coarse, points_test$landuse_predicted_gbm)

## Use GBM for predicted land use
points_test$landuse_predicted <- factor(points_test$landuse_predicted_gbm)

## Examine transition probabilities in validation_landuse_coarse versus predicted land use
dtable <- data.table(points_test)
setkey(dtable, point_id, year)
dtable[, validation_landuse_coarse_next := c(tail(as.character(validation_landuse_coarse), .N-1), as.character(NA)), by="point_id"]
dtable[, landuse_predicted_next := c(tail(as.character(landuse_predicted), .N-1), as.character(NA)), by="point_id"]  # TODO Factor
dtable[, landuse_predicted_two_periods_ahead := c(tail(as.character(landuse_predicted), .N-2), rep(as.character(NA), 2)), by="point_id"]
dtable[, landuse_predicted_two_periods_ahead := factor(landuse_predicted_two_periods_ahead, levels=levels(landuse_predicted))]
head(subset(dtable, select=c("year", "point_id", "validation_landuse_coarse", "validation_landuse_coarse_next")), 30)  # Sanity check

pr_transition <- with(dtable, prop.table(table(validation_landuse_coarse, validation_landuse_coarse_next), 1))
pr_transition_predictions <- with(dtable, prop.table(table(landuse_predicted, landuse_predicted_next), 1))
message("pr_transition: ", paste(round(pr_transition, 3), collapse=" "))
message("pr_transition_predictions: ", paste(round(pr_transition_predictions, 3), collapse=" "))

list_of_initial_hmm_params <- get_initial_hmm_params(opt$landuse_set)

panel <- get_hmm_panel_from_points(dtable, discrete_y_varname="landuse_predicted")  # List of panel elements
baum_welch_time_homogeneous(panel[[1]], params=list_of_initial_hmm_params[[1]])  # Test, compare pi to panel[[1]]$y and panel[[1]]$validation_landuse

list_of_hmm_params_hat <- lapply(list_of_initial_hmm_params, function(initial_hmm_params) {
    return(em_parameter_estimates_time_homogeneous(panel=panel, params=initial_hmm_params, max_iter=50, epsilon=0.001))
})

hmm_params_hat <- list_of_hmm_params_hat[[which.max(sapply(list_of_hmm_params_hat, function(x) max(x$loglik)))]]  # Choose estimate with highest loglik
hmm_params_hat$P  # Compare to pr_transition, pr_transition_predictions
hmm_params_hat$pr_y  # Compare to with(points_test, prop.table(table(validation_landuse_coarse, landuse_predicted), margin=1))
message("hmm_params_hat$P: ", paste(round(hmm_params_hat$P, 3), collapse=" "))

## Run Viterbi and see whether test performance beats raw GBM
viterbi_paths <- lapply(panel, viterbi_path_time_homogeneous, params=hmm_params_hat)
dtable$viterbi_landuse <- levels(dtable$landuse_predicted)[c(viterbi_paths, recursive=TRUE)]
viterbi_test_confusion <- table(dtable$validation_landuse_coarse, dtable$viterbi_landuse)  # Compare to gbm_test_confusion

## Bootstrap panel, compute pr_transition, pr_transition_predictions and HMM estimates on each bootstrap sample
run_bootstrap <- function() {
    panel_indices <- seq_along(panel)
    resampled_panel_indices <- sort(sample(panel_indices, size=length(panel), replace=TRUE))  # Sample by point_id
    panel_boot <- panel[resampled_panel_indices]
    dtable_boot <- rbindlist(lapply(seq_along(panel_boot), function(panel_boot_index) {
        ## Reconstruct dtable from resampled panel; careful, point_id is no longer a unique identifier, use boot_index instead
        panel_element <- panel_boot[[panel_boot_index]]
        data.table(boot_index=panel_boot_index, point_id=panel_element$point_id, y=panel_element$y,
                   validation_landuse_coarse=panel_element$validation_landuse_coarse,
                   validation_landuse=panel_element$validation_landuse, year=seq(min(dtable$year), max(dtable$year)))
    }))
    setkey(dtable_boot, boot_index, year)
    dtable_boot[, predicted_landuse := levels(dtable$landuse_predicted)[y]]
    dtable_boot[, predicted_landuse_next := c(tail(as.character(predicted_landuse), .N-1), as.character(NA)), by="boot_index"]
    dtable_boot[, validation_landuse_coarse_next := c(tail(as.character(validation_landuse_coarse), .N-1), as.character(NA)), by="boot_index"]
    pr_transition_boot <- with(dtable_boot, prop.table(table(validation_landuse_coarse, validation_landuse_coarse_next), 1))
    pr_transition_boot_predictions <- with(dtable_boot, prop.table(table(predicted_landuse, predicted_landuse_next), 1))
    prediction_confusion_matrix <- with(dtable_boot, prop.table(table(validation_landuse_coarse, predicted_landuse), margin=1))  # Compare to hmm_params_hat$pr_y
    list_of_hmm_params_hat <- lapply(list_of_initial_hmm_params, function(initial_hmm_params) {
        return(em_parameter_estimates_time_homogeneous(panel=panel_boot, params=initial_hmm_params, max_iter=50, epsilon=0.001))
    })
    hmm_params_hat <- list_of_hmm_params_hat[[which.max(sapply(list_of_hmm_params_hat, function(x) max(x$loglik)))]]  # Choose estimate with highest loglik
    hmm_params_hat$landuses <- levels(dtable$landuse_predicted)
    return(list(pr_transition=pr_transition_boot,
                pr_transition_predictions=pr_transition_boot_predictions,
                hmm_params_hat=hmm_params_hat,
                prediction_confusion_matrix=prediction_confusion_matrix))
}

n_boostrap_samples <- 100
boots_filename <- sprintf("validation_bootstrap_%s_panel_%s_landuse_set_%s_replications_%s.rds",
                          opt$panel_length, opt$landuse_set, n_boostrap_samples, digest(points_train, algo="crc32"))

boots <- replicate(n_boostrap_samples, run_bootstrap(), simplify=FALSE)
saveRDS(boots, file=boots_filename)

df_boots <- data.frame(replication_index=seq_along(boots))
landuses <- unique(dtable$landuse_predicted_gbm)
for(landuse in landuses) {
    landuse_underscore <- gsub(" |-", "_", landuse)
    varname_pr_transition <- sprintf("pr_%s_%s", landuse_underscore, landuse_underscore)
    varname_pr_transition_predictions <- sprintf("predictions_pr_%s_%s", landuse_underscore, landuse_underscore)
    varname_pr_transition_hmm <- sprintf("hmm_pr_%s_%s", landuse_underscore, landuse_underscore)
    df_boots[, varname_pr_transition] <- sapply(boots, function(x) x$pr_transition[landuse, landuse])
    df_boots[, varname_pr_transition_predictions] <- sapply(boots, function(x) x$pr_transition_predictions[landuse, landuse])
    df_boots[, varname_pr_transition_hmm] <- sapply(boots, function(x) x$hmm_params$P[which(x$hmm_params$landuses == landuse),
                                                                                      which(x$hmm_params$landuses == landuse)])
    test_error <- sprintf("test_error_%s", landuse_underscore)
    hmm_error  <- sprintf("hmm_misclassification_%s", landuse_underscore)
    df_boots[, test_error] <- sapply(boots, function(x) {
        return(1 - x$prediction_confusion_matrix[landuse, landuse])
    })
    df_boots[, hmm_error] <- sapply(boots, function(x) {
        return(1 - x$hmm_params$pr_y[which(x$hmm_params$landuses == landuse),
                                     which(x$hmm_params$landuses == landuse)])
    })
    ## Single transition probability plot showing both HMM and GBM
    df_melted <- melt(subset(df_boots, select=c("replication_index", varname_pr_transition, varname_pr_transition_hmm, varname_pr_transition_predictions)),
                      id.vars=c("replication_index", varname_pr_transition))
    df_melted$label <- ifelse(grepl("^hmm", df_melted$variable), "HMM", "GBM predictions")
    p <- (ggplot(df_melted, aes_string(x=varname_pr_transition, y="value")) +
          geom_point() +
          geom_abline(slope=1, intercept=0, linetype=2, color="grey", size=1.2) +
          xlab(sprintf("probability of %s to %s transition in bootstrapped test data", landuse, landuse)) +
          ylab("transition probability estimate") +
          facet_wrap(~ label, ncol=2))
    outfile <- sprintf("validation_%s_panel_%s_bootstrap_%s_transition_hmm_and_gbm_predictions.png",
                       opt$panel_length, opt$landuse_set, landuse_underscore, landuse_underscore)
    ggsave(outfile, p, width=10, height=8)
    ## GBM test errors and HMM pr_y
    p <- (ggplot(df_boots, aes_string(x=test_error, y=hmm_error)) +
          geom_abline(slope=1, intercept=0, linetype=2, color="grey", size=1.2) +
          xlab(sprintf("error rate for land use predictions in bootstrapped test data, given true land use is %s", landuse)) +
          ylab("HMM estimate of error rate") +
          geom_point())
    outfile <- sprintf("validation_%s_panel_%s_bootstrap_%s_misclassification_probability.png",
                       opt$panel_length, opt$landuse_set, landuse_underscore)
    ggsave(outfile, p, width=10, height=8)
}
outfile <- sprintf("validation_bootstrap_%s_panel_%s_landuse_set_%s_replications.csv",
                   opt$panel_length, opt$landuse_set, n_boostrap_samples)
message("saving ", outfile)
write.csv(df_boots, file=outfile, row.names=FALSE)

message("mean ground truth Pr[S_it+1 = crops | S_it = crops] in boots: ", mean(df_boots$pr_crops_crops))
message("mean ground truth Pr[S_it+1 = pasture | S_it = pasture] in boots: ", mean(df_boots$pr_pasture_pasture))

message("sd ground truth Pr[S_it+1 = crops | S_it = crops] in boots: ", sd(df_boots$pr_crops_crops))
message("sd ground truth Pr[S_it+1 = pasture | S_it = pasture] in boots: ", sd(df_boots$pr_pasture_pasture))

message("mean HMM Pr[S_it+1 = crops | S_it = crops] in boots: ", mean(df_boots$hmm_pr_crops_crops))
message("mean HMM Pr[S_it+1 = pasture | S_it = pasture] in boots: ", mean(df_boots$hmm_pr_pasture_pasture))

message("sd HMM Pr[S_it+1 = crops | S_it = crops] in boots: ", sd(df_boots$hmm_pr_crops_crops))
message("sd HMM Pr[S_it+1 = pasture | S_it = pasture] in boots: ", sd(df_boots$hmm_pr_pasture_pasture))

message("mean Y_it frequency estimate Pr[S_it+1 = crops | S_it = crops] in boots: ", mean(df_boots$predictions_pr_crops_crops))
message("mean Y_it frequency estimate Pr[S_it+1 = pasture | S_it = pasture] in boots: ", mean(df_boots$predictions_pr_pasture_pasture))

message("sd Y_it frequency estimate Pr[S_it+1 = crops | S_it = crops] in boots: ", sd(df_boots$predictions_pr_crops_crops))
message("sd Y_it frequency estimate Pr[S_it+1 = pasture | S_it = pasture] in boots: ", sd(df_boots$predictions_pr_pasture_pasture))

rmse_frequency_pr_pasture_pasture <- with(df_boots, sqrt(mean((predictions_pr_pasture_pasture - pr_pasture_pasture) ^ 2)))
rmse_hmm_pr_pasture_pasture <- with(df_boots, sqrt(mean((hmm_pr_pasture_pasture - pr_pasture_pasture) ^ 2)))
message("RMSEs for Pr[S_it+1 = pasture | S_it = pasture]: HMM ", rmse_hmm_pr_pasture_pasture, ", frequency estimator ", rmse_frequency_pr_pasture_pasture)

rmse_frequency_pr_crops_crops <- with(df_boots, sqrt(mean((predictions_pr_crops_crops - pr_crops_crops) ^ 2)))
rmse_hmm_pr_crops_crops <- with(df_boots, sqrt(mean((hmm_pr_crops_crops - pr_crops_crops) ^ 2)))
message("RMSEs for Pr[S_it+1 = crops | S_it = crops]: HMM ", rmse_hmm_pr_crops_crops, ", frequency estimator ", rmse_frequency_pr_crops_crops)

rmse_hmm_miclassification_crops <- with(df_boots, sqrt(mean((hmm_misclassification_crops - test_error_crops) ^ 2)))
rmse_hmm_miclassification_pasture <- with(df_boots, sqrt(mean((hmm_misclassification_pasture - test_error_pasture) ^ 2)))
message("RMSEs for HMM estimates of misclassification probabilities (i.e. Pr[Y_it | S_it]): crops ", rmse_hmm_miclassification_crops, ", pasture ", rmse_hmm_miclassification_pasture)

message("GBM test set confusion matrix:")
gbm_test_confusion
