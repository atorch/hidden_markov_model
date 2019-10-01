library(data.table)
library(gbm)
library(ggplot2)
library(optparse)
library(plyr); library(dplyr)  # For bind_rows
library(randomForest)
library(sp)

set.seed(123123)

opt_list <- list(make_option("--landuse_set", default="binary"))  # Binary for {crops, pasture}, ternary for {soy, other crops, pasture}
opt <- parse_args(OptionParser(option_list=opt_list))
message("command line options: ", paste(sprintf("%s=%s", names(opt), opt), collapse=", "))
stopifnot(opt$landuse_set %in% c("binary", "ternary"))

source("embrapa_validation_landuse_mapping.R")
source("hmm_functions.R")
source("ggplot_utils.R")
set_ggplot_theme()

points <- readRDS("~/Dropbox/amazon_hmm_shared/embrapa_validation/embrapa_validation_points_with_covariates.rds")
message("done loading ", nrow(points), " point-years, ", length(unique(points$point_id)), " unique spatial points")

points_with_validation_landuse <- readRDS("~/Dropbox/amazon_hmm_shared/embrapa_validation/embrapa_validation_points_with_validation_landuse.rds")
points <- merge(points, points_with_validation_landuse, all.x=TRUE, by=c("point_id", "year"))  # Adds validation_landuse column
points <- arrange(points, point_id, year)

state_abbr <- unique(as.character(points$state))  # E.g. MT for Mato Grosso
stopifnot(length(state_abbr) == 1)
stopifnot(nchar(state_abbr) == 2)
message("all points are in ", state_abbr)

validation_years <- sort(unique(points$year[!is.na(points$validation_landuse)]))
message("years with validation land use data: ", paste(validation_years, collapse=", "))

## Short panel: keep only validation years (i.e. drop 2001-2005 and 2011-2016)
points <- subset(points, year %in% validation_years)
message("short panel, keep only ", nrow(points), " point-years, ", length(unique(points$point_id)), " unique spatial points")

if(opt$landuse_set == "binary") {
    points$validation_landuse_coarse <- plyr::revalue(points$validation_landuse, replace=map_landuse_to_S_binary)
    stopifnot(all(is.na(points$validation_landuse_coarse) | points$validation_landuse_coarse %in% c("crops", "pasture")))
} else {
    points$validation_landuse_coarse <- plyr::revalue(points$validation_landuse, replace=map_landuse_to_S_ternary)
    stopifnot(all(is.na(points$validation_landuse_coarse) | points$validation_landuse_coarse %in% c("soy", "non-soy crops", "pasture")))
}

table(points$validation_landuse_coarse)  # Careful, pasture is rare
points$validation_landuse_coarse <- factor(points$validation_landuse_coarse)  # Factor for randomForest classification

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

point_id_with_validation_landuse_coarse <- unique(points$point_id[!is.na(points$validation_landuse_coarse)])
points <- subset(points, point_id %in% point_id_with_validation_landuse_coarse)  # Drop points whose validation landuse is always NA

point_id_train <- sample(unique(points$point_id), size=floor(0.60 * length(unique(points$point_id))))
points_train <- subset(points, point_id %in% point_id_train)
points_test <- subset(points, !point_id %in% point_id_train)
table(points_train$validation_landuse_coarse)  # Careful, pasture is rare -- should make sure that training set is not 100% crops, same for test
table(points_test$validation_landuse_coarse)  # Similar comment about non-soy crops in ternary case

predictors <- c(evi_vars, nir_vars, mir_vars, red_vars, blue_vars)
model_formula <- reformulate(termlabels=predictors, response="validation_landuse_coarse")
model <- randomForest(formula=model_formula, ntree=1000, nodesize=1, do.trace=TRUE,
                      mtry=floor(sqrt(length(predictors))), na.action=na.omit, data=points_train)
model$confusion  # Large OOB error rate for pasture -- non-soy crops even worse in ternary case

points_test$landuse_predicted_rf <- predict(model, type="response", newdata=points_test)
mean(is.na(points_test$landuse_predicted_rf))  # Around 0.11 -- missing when MODIS variables are NA
table(points_test$landuse_predicted_rf)
table(points_test$validation_landuse_coarse, points_test$landuse_predicted_rf)  # RF test confusion matrix, compare to OOB

## Fit GBM, see whether it beats RF -- takes around 10 minutes to fit -- appears to beat RF for rare landuse in both binary and ternary case
model_gbm <- gbm(formula=model_formula, distribution="multinomial",  # Does this have to be bernoulli for two-class case?
                 data=subset(points, !is.na(validation_landuse_coarse)), n.trees=8000, interaction.depth=3, shrinkage=0.001, cv.folds=3)
ntrees_star <- gbm.perf(model_gbm, method="cv")  # Plot shows training and CV deviance as function of number of trees -- selects around 5000-6000 trees
gbm_predictions <- predict(model_gbm, type="response", newdata=points_test, ntrees=ntrees_star)  # Outputs class probabilities
points_test$landuse_predicted_gbm <- colnames(gbm_predictions)[as.integer(apply(gbm_predictions, 1, which.max))]
table(points_test$landuse_predicted_gbm)
table(points_test$validation_landuse_coarse, points_test$landuse_predicted_gbm)  # GBM test confusion matrix, compare to RF
table(points_test$landuse_predicted_gbm, points_test$landuse_predicted_rf)  # GBM versus RF predictions

## Use GBM for predicted land use
points_test$landuse_predicted <- factor(points_test$landuse_predicted_gbm)

## Examine transition probabilities in validation_landuse_coarse versus predicted land use
dtable <- data.table(points_test)
setkey(dtable, point_id, year)
dtable[, validation_landuse_coarse_next := c(tail(as.character(validation_landuse_coarse), .N-1), as.character(NA)), by="point_id"]
dtable[, landuse_predicted_next := c(tail(as.character(landuse_predicted), .N-1), as.character(NA)), by="point_id"]
head(subset(dtable, select=c("year", "point_id", "validation_landuse_coarse", "validation_landuse_coarse_next")), 30)  # Sanity check

pr_transition <- with(dtable, prop.table(table(validation_landuse_coarse, validation_landuse_coarse_next), 1))
pr_transition_predictions <- with(dtable, prop.table(table(landuse_predicted, landuse_predicted_next), 1))

if(opt$landuse_set == "binary") {
    list_of_initial_hmm_params <- list(list(P=rbind(c(0.90, 0.10),
                                                    c(0.10, 0.90)),
                                            mu=c(0.5, 0.5),
                                            pr_y=rbind(c(0.90, 0.10),
                                                       c(0.10, 0.90)),
                                            n_components=2),
                                       list(P=rbind(c(0.80, 0.20),
                                                    c(0.20, 0.80)),
                                            mu=c(0.6, 0.4),
                                            pr_y=rbind(c(0.90, 0.10),
                                                       c(0.20, 0.80)),
                                            n_components=2),
                                       list(P=rbind(c(0.90, 0.10),
                                                    c(0.50, 0.50)),
                                            mu=c(0.5, 0.5),
                                            pr_y=rbind(c(0.90, 0.10),
                                                       c(0.50, 0.50)),
                                            n_components=2))
} else {
    list_of_initial_hmm_params <- list(list(P=rbind(c(0.80, 0.10, 0.10),
                                                    c(0.10, 0.80, 0.10),
                                                    c(0.10, 0.10, 0.80)),
                                            mu=c(1/3, 1/3, 1/3),
                                            pr_y=rbind(c(0.80, 0.10, 0.10),
                                                       c(0.10, 0.80, 0.10),
                                                       c(0.10, 0.10, 0.80)),
                                            n_components=3),
                                       list(P=rbind(c(.50, 0.25, 0.25),
                                                    c(0.25, 0.50, 0.25),
                                                    c(0.10, 0.10, 0.80)),
                                            mu=c(0.2, 0.3, 0.5),
                                            pr_y=rbind(c(0.40, 0.20, 0.40),
                                                       c(0.20, 0.40, 0.40),
                                                       c(0.10, 0.10, 0.80)),
                                            n_components=3))
}

panel <- get_hmm_panel_from_points(dtable, discrete_y_varname="landuse_predicted")  # List of panel elements
baum_welch_time_homogeneous(panel[[1]], params=list_of_initial_hmm_params[[1]])  # Test, compare pi to panel[[1]]$y and panel[[1]]$validation_landuse

list_of_hmm_params_hat <- lapply(list_of_initial_hmm_params, function(initial_hmm_params) {
    return(em_parameter_estimates_time_homogeneous(panel=panel, params=initial_hmm_params, max_iter=50, epsilon=0.001))
})

hmm_params_hat <- list_of_hmm_params_hat[[which.max(sapply(list_of_hmm_params_hat, function(x) max(x$loglik)))]]  # Choose estimate with highest loglik
hmm_params_hat$P  # Compre to pr_transition, pr_transition_predictions
hmm_params_hat$pr_y  # Compare to with(points_test, prop.table(table(validation_landuse_coarse, landuse_predicted), margin=1))

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
boots <- replicate(n_boostrap_samples, run_bootstrap(), simplify=FALSE)
outfile <- sprintf("bootstrap_%s_landuse_set_%s_replications.rds", opt$landuse_set, n_boostrap_samples)
saveRDS(boots, file=outfile)

## TODO Lines below assume we're in the binary {crops, pasture} case, need to edit to handle opt$landuse_set == "ternary"
df_boots <- data.frame(replication_index=seq_along(boots),
                       pr_crops_crops=sapply(boots, function(x) x$pr_transition["crops", "crops"]),
                       pr_pasture_pasture=sapply(boots, function(x) x$pr_transition["pasture", "pasture"]),
                       predictions_pr_crops_crops=sapply(boots, function(x) x$pr_transition_predictions["crops", "crops"]),
                       predictions_pr_pasture_pasture=sapply(boots, function(x) x$pr_transition_predictions["pasture", "pasture"]),
                       hmm_pr_crops_crops=sapply(boots, function(x) x$hmm_params$P[which(x$hmm_params$landuses == "crops"),
                                                                                   which(x$hmm_params$landuses == "crops")]),
                       hmm_pr_pasture_pasture=sapply(boots, function(x) x$hmm_params$P[which(x$hmm_params$landuses == "pasture"),
                                                                                       which(x$hmm_params$landuses == "pasture")]),
                       test_error_crops=sapply(boots, function(x) x$prediction_confusion_matrix["crops", "pasture"]),
                       test_error_pasture=sapply(boots, function(x) x$prediction_confusion_matrix["pasture", "crops"]),
                       hmm_pr_y_crops=sapply(boots, function(x) x$hmm_params$pr_y[which(x$hmm_params$landuses == "crops"),
                                                                                  which(x$hmm_params$landuses == "pasture")]),
                       hmm_pr_y_pasture=sapply(boots, function(x) x$hmm_params$pr_y[which(x$hmm_params$landuses == "pasture"),
                                                                                    which(x$hmm_params$landuses == "crops")]))
outfile <- sprintf("~/Dropbox/amazon_hmm_shared/embrapa_validation/bootstrap_%s_landuse_set_%s_replications.csv", opt$landuse_set, n_boostrap_samples)
write.csv(df_boots, file=outfile, row.names=FALSE)

p <- (ggplot(df_boots, aes(x=pr_pasture_pasture, y=hmm_pr_pasture_pasture)) +
      geom_abline(slope=1, intercept=0, linetype=2, color="grey", size=1.2) +
      xlab("probability of pasture-to-pasture transition in bootstrapped test data") +
      ylab("HMM transition probability estimate") +
      geom_point())
p
ggsave("~/Dropbox/amazon_hmm_shared/plots/embrapa_validation_bootstrap_pasture_transition_hmm.png", p, width=10, height=8)

p <- (ggplot(df_boots, aes(x=pr_pasture_pasture, y=predictions_pr_pasture_pasture)) +
      geom_abline(slope=1, intercept=0, linetype=2, color="grey", size=1.2) +
      xlab("probability of pasture-to-pasture transition in bootstrapped test data") +
      ylab("transition probability estimate from land use predictions") +
      geom_point())
p
ggsave("~/Dropbox/amazon_hmm_shared/plots/embrapa_validation_bootstrap_pasture_transition_predictions.png", p, width=10, height=8)

p <- (ggplot(df_boots, aes(x=pr_crops_crops, y=hmm_pr_crops_crops)) +
      geom_abline(slope=1, intercept=0, linetype=2, color="grey", size=1.2) +
      xlab("probability of crops-to-crops transition in bootstrapped test data") +
      ylab("HMM transition probability estimate") +
      geom_point())
p
ggsave("~/Dropbox/amazon_hmm_shared/plots/embrapa_validation_bootstrap_crops_transition_hmm.png", p, width=10, height=8)

p <- (ggplot(df_boots, aes(x=pr_crops_crops, y=predictions_pr_crops_crops)) +
      geom_abline(slope=1, intercept=0, linetype=2, color="grey", size=1.2) +
      xlab("probability of crops-to-crops transition in bootstrapped test data") +
      ylab("transition probability estimate from land use predictions") +
      geom_point())
p
ggsave("~/Dropbox/amazon_hmm_shared/plots/embrapa_validation_bootstrap_crops_transition_predictions.png", p, width=10, height=8)

p <- (ggplot(df_boots, aes(x=test_error_crops, y=hmm_pr_y_crops)) +
      geom_abline(slope=1, intercept=0, linetype=2, color="grey", size=1.2) +
      xlab("error rate for land use predictions in bootstrapped test data, given true land use is crops") +
      ylab("HMM estimate of error rate") +
      geom_point())
p
ggsave("~/Dropbox/amazon_hmm_shared/plots/embrapa_validation_bootstrap_crops_misclassification_probability.png", p, width=10, height=8)

p <- (ggplot(df_boots, aes(x=test_error_pasture, y=hmm_pr_y_pasture)) +
      geom_abline(slope=1, intercept=0, linetype=2, color="grey", size=1.2) +
      xlab("error rate for land use predictions in bootstrapped test data, given true land use is pasture") +
      ylab("HMM estimate of error rate") +
      geom_point())
p  # HMM estimate of pasture error appears to be biased downward
ggsave("~/Dropbox/amazon_hmm_shared/plots/embrapa_validation_bootstrap_pasture_misclassification_probability.png", p, width=10, height=8)
