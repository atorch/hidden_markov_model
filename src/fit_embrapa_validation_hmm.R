library(data.table)
library(digest)
library(gbm)
library(ggplot2)
library(optparse)
library(parallel)
library(plyr); library(dplyr)  # For bind_rows
library(randomForest)
library(Rsolnp)
library(sp)
library(stringr)

set.seed(789)

opt_list <- list(make_option("--n_bootstrap_samples", default=100, type=("integer")),
                 make_option("--panel_length", default="short"),
                 make_option("--classifier_training_fraction", default=0.15, type="double"),
                 make_option("--classifier_pasture_fraction", default=0.6, type="double",
                             help="This parameter is used to oversample pasture (to make the GBM training set more balanced)."))
opt <- parse_args(OptionParser(option_list=opt_list))
message("command line options: ", paste(sprintf("%s=%s", names(opt), opt), collapse=", "))

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
    ## Estimate HMM using "long" panel, i.e. including years for which we don't observe ground truth land cover
    points <- subset(points, year >= min(validation_years))
    
} else {
    ## Estimate HMM using only those years for which we observe ground truth land cover
    points <- subset(points, year %in% validation_years)
}

message("keeping the following years: ", paste(unique(points$year), collapse=", "))
message("keeping ", nrow(points), " point-years, ", length(unique(points$point_id)), " unique spatial points")

point_id_ever_mata <- unique(points$point_id[!is.na(points$validation_landuse) & tolower(points$validation_landuse) == "mata"])
point_id_ever_reflorestamento <- unique(points$point_id[!is.na(points$validation_landuse) & tolower(points$validation_landuse) == "reflorestamento"])
points <- subset(points, !point_id %in% c(point_id_ever_mata, point_id_ever_reflorestamento))  # Down to 403 points
message("keeping ", nrow(points), " point-years, ", length(unique(points$point_id)), " unique spatial points after removing points whose validation landuse was ever mata or reflorestamento")

message("number of point-years missing ground truth land use: ", sum(is.na(points$validation_landuse)))
message("fraction of point-years missing ground truth land use: ", mean(is.na(points$validation_landuse)))

message("number of unique points missing ground truth land use in one or more years: ", length(unique(points$point_id[is.na(points$validation_landuse)])))

## Binary for {crops, pasture}, ternary for {soy, other crops, pasture}
map_landuse_to_S <- map_landuse_to_S_binary
points$validation_landuse_coarse <- plyr::revalue(points$validation_landuse, replace=map_landuse_to_S)
stopifnot(all(is.na(points$validation_landuse_coarse) | points$validation_landuse_coarse %in% c("crops", "pasture")))

table(points$validation_landuse_coarse)  # Careful, pasture is rare
points$validation_landuse_coarse <- factor(points$validation_landuse_coarse)  # Factor for classification

## Note: pasture would be a small percentage of the training set if we took a simple random sample,
## so we stratify in order to over-sample pasture for the GBM training set
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

message("GBM test set accuracy: ", mean(points_test$landuse_predicted_gbm == points_test$validation_landuse_coarse, na.rm=TRUE))

## Use GBM for predicted land use
points_test$landuse_predicted <- factor(points_test$landuse_predicted_gbm)
pr_Y_given_S <- with(points_test, prop.table(table(validation_landuse_coarse, landuse_predicted), margin=1))

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

## Note: time is needed for time-varying estimation code
dtable[, time := year - min(year)]
panel <- get_hmm_panel_from_points(dtable, discrete_y_varname="landuse_predicted")  # List of panel elements

list_of_hmm_params_hat <- lapply(initial_hmm_params, function(initial_hmm_params) {
    return(em_parameter_estimates_time_homogeneous(panel=panel, params=initial_hmm_params, max_iter=50, epsilon=0.001))
})

## Choose estimate with highest log likelihood
hmm_params_hat <- list_of_hmm_params_hat[[which.max(sapply(list_of_hmm_params_hat, function(x) max(x$loglik)))]]

hmm_params_hat$P  # Compare to pr_transition, pr_transition_predictions
hmm_params_hat$pr_y  # Compare to pr_Y_given_S
message("hmm_params_hat$P: ", paste(round(hmm_params_hat$P, 3), collapse=" "))

## Run Viterbi and see whether test performance beats raw GBM
viterbi_paths <- lapply(panel, viterbi_path_time_homogeneous, params=hmm_params_hat)
dtable$viterbi_landuse <- levels(dtable$landuse_predicted)[c(viterbi_paths, recursive=TRUE)]
viterbi_test_confusion <- table(dtable$validation_landuse_coarse, dtable$viterbi_landuse)  # Compare to gbm_test_confusion

## TODO Include (viterbi_test_confusion - gbm_test_confusion) in the validation section?
## HMM improves classification accuracy "for free"! We get more mass on the diagonal of the confusion matrix
## TODO Include improvement in accuracy from using viterbi in boostrap?

## Run MD with time-homogeneous parameters
## TODO Make this function also return EM estimates, move code into hmm_functions.R
estimates_time_homogeneous <- get_minimum_distance_estimates_random_initialization_time_homogeneous(initial_hmm_params[[1]], panel)

estimates_time_homogeneous$min_dist_params_hat_best_objfn

## Try running both MD and EM with time-varying parameters
params0 <- initial_hmm_params[[1]]
params0$P_list <- rep(list(params0$P), length(unique((dtable$year))) - 1)
params0$P <- NULL

estimates_time_varying <- get_hmm_and_minimum_distance_estimates_random_initialization(params0, panel)

## Estimated crops->pasture transition probability hits edge (it's zero) in one period, but seems reasonable
estimates_time_varying$min_dist_params_hat_best_objfn

## Also looks good!  Doesn't hit edge of parameter space
estimates_time_varying$em_params_hat_best_likelihood

## Compare time-varying MD and EM estimates (they're fairly close: mu and pr_y are nearly identical, P_list differs in some periods)
estimates_time_varying$em_params_hat_best_likelihood$pr_y - estimates_time_varying$min_dist_params_hat_best_objfn$pr_y
estimates_time_varying$em_params_hat_best_likelihood$mu - estimates_time_varying$min_dist_params_hat_best_objfn$mu
lapply(seq_along(params0$P_list), function(time) {
    estimates_time_varying$em_params_hat_best_likelihood$P_list[[time]] - estimates_time_varying$min_dist_params_hat_best_objfn$P_list[[time]]
})

## Bootstrap panel, compute pr_transition, pr_transition_predictions and HMM estimates on each bootstrap sample
run_bootstrap <- function() {
    require(data.table)
    panel_indices <- seq_along(panel)
    resampled_panel_indices <- sort(sample(panel_indices, size=length(panel), replace=TRUE))  # Sample by point_id

    panel_boot <- panel[resampled_panel_indices]
    dtable_boot <- rbindlist(lapply(seq_along(panel_boot), function(panel_boot_index) {
        ## Reconstruct dtable from resampled panel
        panel_element <- panel_boot[[panel_boot_index]]
        data.table(boot_index=panel_boot_index, point_id=panel_element$point_id, y=panel_element$y,
                   validation_landuse_coarse=panel_element$validation_landuse_coarse,
                   validation_landuse=panel_element$validation_landuse, year=seq(min(dtable$year), max(dtable$year)))
    }))

    ## Careful, point_id is no longer a unique identifier (because of sampling with replacement), so we use boot_index instead
    setkey(dtable_boot, boot_index, year)

    dtable_boot[, predicted_landuse := levels(dtable$landuse_predicted)[y]]
    dtable_boot[, predicted_landuse_next := c(tail(as.character(predicted_landuse), .N-1), as.character(NA)), by="boot_index"]
    dtable_boot[, validation_landuse_coarse_next := c(tail(as.character(validation_landuse_coarse), .N-1), as.character(NA)), by="boot_index"]

    pr_transition_boot <- with(dtable_boot, prop.table(table(validation_landuse_coarse, validation_landuse_coarse_next), 1))
    pr_transition_boot_predictions <- with(dtable_boot, prop.table(table(predicted_landuse, predicted_landuse_next), 1))

    prediction_confusion_matrix <- with(dtable_boot, prop.table(table(validation_landuse_coarse, predicted_landuse), margin=1))

    list_of_hmm_params_hat <- lapply(initial_hmm_params, function(initial_hmm_params) {
        return(em_parameter_estimates_time_homogeneous(panel=panel_boot, params=initial_hmm_params, max_iter=50, epsilon=0.001))
    })

    ## Choose estimate with highest loglik
    hmm_params_hat <- list_of_hmm_params_hat[[which.max(sapply(list_of_hmm_params_hat, function(x) max(x$loglik)))]]

    estimates_time_homogeneous <- get_minimum_distance_estimates_random_initialization_time_homogeneous(initial_hmm_params[[1]], panel_boot)

    md_params_hat <- estimates_time_homogeneous$min_dist_params_hat_best_objfn

    params0 <- initial_hmm_params[[1]]
    params0$P_list <- rep(list(params0$P), length(unique((dtable$year))) - 1)
    params0$P <- NULL

    estimates_time_varying <- get_hmm_and_minimum_distance_estimates_random_initialization(params0, panel_boot)

    hmm_params_hat_time_varying <- estimates_time_varying$em_params_hat_best_likelihood
    md_params_hat_time_varying <- estimates_time_varying$min_dist_params_hat_best_objfn

    ## Run Viterbi and see whether test performance beats raw GBM
    ## TODO Both EM and MD
    viterbi_paths_time_varying <- lapply(panel_boot, viterbi_path, params=hmm_params_hat_time_varying)
    dtable_boot$viterbi_landuse_time_varying <- levels(dtable$landuse_predicted)[c(viterbi_paths_time_varying, recursive=TRUE)]

    viterbi_paths <- lapply(panel_boot, viterbi_path_time_homogeneous, params=hmm_params_hat)
    dtable_boot$viterbi_landuse <- levels(dtable$landuse_predicted)[c(viterbi_paths, recursive=TRUE)]

    gbm_accuracy <- mean(dtable_boot$predicted_landuse == dtable_boot$validation_landuse_coarse, na.rm=TRUE)
    viterbi_accuracy_time_varying <- mean(dtable_boot$viterbi_landuse_time_varying == dtable_boot$validation_landuse_coarse, na.rm=TRUE)
    viterbi_accuracy <- mean(dtable_boot$viterbi_landuse == dtable_boot$validation_landuse_coarse, na.rm=TRUE)

    ## TODO May need to adjust this dataframe to make the Classification Accuracy plot easier to read
    df_classification_accuracy <- data.frame("variable"="Classification Accuracy",
                                             "transition_year"="All Years",
                                             "estimator_type"=c("", "Time Homogeneous", "Time Varying"),
                                             "estimator"=c("GBM",
                                                           "GBM with Time Homogeneous HMM Correction (Viterbi)",
                                                           "GBM with Time Varying HMM Correction (Viterbi)"),
                                             "estimated_value"=c(gbm_accuracy, viterbi_accuracy, viterbi_accuracy_time_varying))

    transition_years <- seq(min(dtable$year), max(dtable$year) - 1)

    pr_transition_boot_time_varying <- lapply(transition_years, function(fixed_year) {
        with(subset(dtable_boot, year == fixed_year), prop.table(table(validation_landuse_coarse, validation_landuse_coarse_next), 1))
    })
    pr_transition_boot_predictions_time_varying <- lapply(transition_years, function(fixed_year) {
        with(subset(dtable_boot, year == fixed_year), prop.table(table(predicted_landuse, predicted_landuse_next), 1))
    })

    ## TODO Run viterbi with time-varying transition probs, record improvement in classification accuracy

    df_pr_crops_pasture <- data.frame("variable"="Crops to Pasture",
                                      "transition_year"="All Years",
                                      "estimator_type"="Time Homogeneous",
                                      "estimator"=c("EM", "MD", "Frequency", "Ground Truth"),
                                      "estimated_value"=c(hmm_params_hat$P[1, 2],
                                                          md_params_hat$P[1, 2],
                                                          pr_transition_boot_predictions[1, 2],
                                                          pr_transition_boot[1, 2]))

    df_pr_pasture_crops <- data.frame("variable"="Pasture to Crops",
                                      "transition_year"="All Years",
                                      "estimator_type"="Time Homogeneous",
                                      "estimator"=c("EM", "MD", "Frequency", "Ground Truth"),
                                      "estimated_value"=c(hmm_params_hat$P[2, 1],
                                                          md_params_hat$P[2, 1],
                                                          pr_transition_boot_predictions[2, 1],
                                                          pr_transition_boot[2, 1]))

    df_pr_pasture_crops_time_varying_errors <- do.call(rbind, lapply(transition_years, function(fixed_year) {
        index <- which(transition_years == fixed_year)
        data.frame("variable"="Pasture to Crops",
                   "transition_year"=paste0(fixed_year, "-", fixed_year + 1),
                   "estimator_type"="Time Varying (Errors Relative to Ground Truth)",
                   "estimator"=c("EM", "MD", "Frequency"),
                   "estimated_value"=c(hmm_params_hat_time_varying$P_list[[index]][2, 1] - pr_transition_boot_time_varying[[index]][2, 1],
                                       md_params_hat_time_varying$P_list[[index]][2, 1] - pr_transition_boot_time_varying[[index]][2, 1],
                                       pr_transition_boot_predictions_time_varying[[index]][2, 1] - pr_transition_boot_time_varying[[index]][2, 1]))
    }))

    df_pr_pasture_crops_time_varying <- do.call(rbind, lapply(transition_years, function(fixed_year) {
        index <- which(transition_years == fixed_year)
        data.frame("variable"="Pasture to Crops",
                   "transition_year"=paste0(fixed_year, "-", fixed_year + 1),
                   "estimator_type"="Time Varying",
                   "estimator"=c("EM", "MD", "Frequency", "Ground Truth"),
                   "estimated_value"=c(hmm_params_hat_time_varying$P_list[[index]][2, 1],
                                       md_params_hat_time_varying$P_list[[index]][2, 1],
                                       pr_transition_boot_predictions_time_varying[[index]][2, 1],
                                       pr_transition_boot_time_varying[[index]][2, 1]))
    }))

    df_pr_crops_pasture_time_varying_errors <- do.call(rbind, lapply(transition_years, function(fixed_year) {
        index <- which(transition_years == fixed_year)
        data.frame("variable"="Crops to Pasture",
                   "transition_year"=paste0(fixed_year, "-", fixed_year + 1),
                   "estimator_type"="Time Varying (Errors Relative to Ground Truth)",
                   "estimator"=c("EM", "MD", "Frequency"),
                   "estimated_value"=c(hmm_params_hat_time_varying$P_list[[index]][1, 2] - pr_transition_boot_time_varying[[index]][1, 2],
                                       md_params_hat_time_varying$P_list[[index]][1, 2] - pr_transition_boot_time_varying[[index]][1, 2],
                                       pr_transition_boot_predictions_time_varying[[index]][1, 2] - pr_transition_boot_time_varying[[index]][1, 2]))
    }))

    df_pr_crops_pasture_time_varying <- do.call(rbind, lapply(transition_years, function(fixed_year) {
        index <- which(transition_years == fixed_year)
        data.frame("variable"="Crops to Pasture",
                   "transition_year"=paste0(fixed_year, "-", fixed_year + 1),
                   "estimator_type"="Time Varying",
                   "estimator"=c("EM", "MD", "Frequency", "Ground Truth"),
                   "estimated_value"=c(hmm_params_hat_time_varying$P_list[[index]][1, 2],
                                       md_params_hat_time_varying$P_list[[index]][1, 2],
                                       pr_transition_boot_predictions_time_varying[[index]][1, 2],
                                       pr_transition_boot_time_varying[[index]][1, 2]))
    }))

    ## Estimated Pr[Y=crops | S=crops] for models in which transition probabilities are time-varying
    df_pr_y_crops_time_varying <- data.frame("variable"="Pr[Y = crops | S = crops]",
                                             "transition_year"="",
                                             "estimator_type"="Time Varying",
                                             "estimator"=c("EM", "MD", "Ground Truth"),
                                             "estimated_value"=c(hmm_params_hat_time_varying$pr_y[1, 1],
                                                                 md_params_hat_time_varying$pr_y[1, 1],
                                                                 prediction_confusion_matrix[1, 1]))

    ## Estimated Pr[Y=pasture | S=pasture] for models in which transition probabilities are time-varying
    df_pr_y_pasture_time_varying <- data.frame("variable"="Pr[Y = pasture | S = pasture]",
                                               "transition_year"="",
                                               "estimator_type"="Time Varying",
                                               "estimator"=c("EM", "MD", "Ground Truth"),
                                               "estimated_value"=c(hmm_params_hat_time_varying$pr_y[2, 2],
                                                                   md_params_hat_time_varying$pr_y[2, 2],
                                                                   prediction_confusion_matrix[2, 2]))

    ## TODO Add MD Pr[Y | S]
    df_pr_y_crops <- data.frame("variable"="Pr[Y = crops | S = crops]",
                                "transition_year"="",
                                "estimator_type"="Time Homogeneous",
                                "estimator"=c("EM", "Ground Truth"),
                                "estimated_value"=c(hmm_params_hat$pr_y[1, 1], prediction_confusion_matrix[1, 1]))

    df_pr_y_pasture <- data.frame("variable"="Pr[Y = pasture | S = pasture]",
                                  "transition_year"="",
                                  "estimator_type"="Time Homogeneous",
                                  "estimator"=c("EM", "Ground Truth"),
                                  "estimated_value"=c(hmm_params_hat$pr_y[2, 2], prediction_confusion_matrix[2, 2]))

    return(rbind(df_pr_crops_pasture,
                 df_pr_pasture_crops,
                 df_pr_y_crops,
                 df_pr_y_pasture,
                 df_pr_y_crops_time_varying,
                 df_pr_y_pasture_time_varying,
                 df_pr_pasture_crops_time_varying_errors,
                 df_pr_pasture_crops_time_varying,
                 df_pr_crops_pasture_time_varying_errors,
                 df_pr_crops_pasture_time_varying,
                 df_classification_accuracy))

}

num_cores <- detectCores()
cluster <- makeCluster(num_cores)  # Call stopCluster when done

clusterExport(cluster, c("baum_welch",
                         "baum_welch_time_homogeneous",
                         "dtable",
                         "em_parameter_estimates_time_homogeneous",
                         "eq_function_minimum_distance",
                         "eq_function_min_dist_time_homogeneous",
                         "get_expectation_maximization_estimates",
                         "get_hmm_and_minimum_distance_estimates_random_initialization",
                         "get_min_distance_estimates",
                         "get_min_distance_estimates_time_homogeneous",
                         "get_minimum_distance_estimates_random_initialization_time_homogeneous",
                         "get_random_initial_parameters",
                         "get_transition_probs_from_M_S_joint",
                         "initial_hmm_params",
                         "is_diag_dominant",
                         "objfn_minimum_distance",
                         "objfn_min_dist_time_homogeneous",
                         "panel",
                         "run_bootstrap",
                         "valid_panel_element",
                         "valid_parameters",
                         "valid_parameters_time_homogeneous",
                         "viterbi_path",
                         "viterbi_path_time_homogeneous"))

boots <- parLapply(cluster, seq_len(opt$n_bootstrap_samples), function(unused_input) {
    run_bootstrap()
})
boots <- rbindlist(boots)

stopCluster(cluster)

boots_filename <- sprintf("validation_bootstrap_%s_panel_%s_replications_%s.rds",
                          opt$panel_length, opt$n_bootstrap_samples, digest(points_train, algo="crc32"))
## saveRDS(boots, file=boots_filename)

## Set levels to control order in which estimators appear in plots
boots[, estimator_factor := factor(estimator,
                                   levels=c("Frequency", "EM", "MD", "Ground Truth", "GBM", "GBM with Time Homogeneous HMM Correction (Viterbi)", "GBM with Time Varying HMM Correction (Viterbi)"))]

boots_summary <- boots[, list("mean_estimated_value"=mean(estimated_value),
                              "sd_estimated_value"=sd(estimated_value)),
                       by=c("variable", "estimator", "estimator_factor", "estimator_type", "transition_year")]

## TODO Cutoff at zero for lower bound?  pmax(0, ...)  Some of the CIs for transition probabilities include negative values
boots_summary[, lb := mean_estimated_value - 1.96 * sd_estimated_value]
boots_summary[, ub := mean_estimated_value + 1.96 * sd_estimated_value]

boots_summary_P <- subset(boots_summary, estimator_type == "Time Homogeneous" & variable %in% c("Crops to Pasture", "Pasture to Crops"))
boots_summary_P_time_varying <- subset(boots_summary, estimator_type == "Time Varying" & variable %in% c("Crops to Pasture", "Pasture to Crops"))

boots_summary_P_time_varying_errors <- subset(boots_summary, estimator_type == "Time Varying (Errors Relative to Ground Truth)")

boots_summary_accuracy <- subset(boots_summary, variable == "Classification Accuracy")

p <- ggplot(boots_summary_accuracy,
            aes(x = mean_estimated_value, y = estimator_factor, xmin = lb, xmax = ub)) +
    geom_point() +
    geom_errorbarh(height=0) +
    scale_x_continuous('Classification Accuracy') +
    theme(axis.title.y=element_blank())
ggsave("embrapa_bootstrap_classification_accuracy_confidence_intervals.png", width=10, height=8)

p <- ggplot(boots_summary_P,
            aes(x = mean_estimated_value, y = estimator_factor, xmin = lb, xmax = ub)) +
    geom_point() +
    geom_errorbarh(height=0) +
    scale_x_continuous('Transition Rate') +
    theme(axis.title.y=element_blank()) +
    facet_wrap(~ variable, scales='free_x')
ggsave("embrapa_bootstrap_transition_probability_time_homogeneous_confidence_intervals.png", p, width=10, height=8)

p <- ggplot(boots_summary_P_time_varying_errors,
            aes(x = mean_estimated_value, y = estimator_factor, xmin = lb, xmax = ub)) +
    geom_point() +
    geom_errorbarh(height=0) +
    scale_x_continuous('Error in Estimated Transition Rate') +
    theme(axis.title.y=element_blank()) +
    facet_grid(transition_year ~ variable, scales='free') +
    geom_vline(xintercept = 0, linetype=2)
ggsave("embrapa_bootstrap_transition_probability_time_varying_errors_confidence_intervals.png", p, width=10, height=12)

p <- ggplot(boots_summary_P_time_varying,
            aes(x = mean_estimated_value, y = estimator_factor, xmin = lb, xmax = ub)) +
    geom_point() +
    geom_errorbarh(height=0) +
    scale_x_continuous('Transition Rate') +
    theme(axis.title.y=element_blank()) +
    facet_grid(transition_year ~ variable, scales='free')
ggsave("embrapa_bootstrap_transition_probability_time_varying_confidence_intervals.png", p, width=10, height=12)

boots_summary_pr_y <- subset(boots_summary, variable %in% c("Pr[Y = pasture | S = pasture]", "Pr[Y = crops | S = crops]"))

p <- ggplot(boots_summary_pr_y,
            aes(x = mean_estimated_value, y = estimator_factor, xmin = lb, xmax = ub)) +
    geom_point() +
    geom_errorbarh(height=0) +
    scale_x_continuous('Pr[Y | S]') +
    theme(axis.title.y=element_blank()) +
    facet_grid(estimator_type ~ variable, scales='free_x')
ggsave("embrapa_bootstrap_pr_y_given_s_confidence_intervals.png", p, width=10, height=8)

message("GBM test set confusion matrix:")
print(gbm_test_confusion)

message("GBM test set accuracy: ", mean(points_test$landuse_predicted_gbm == points_test$validation_landuse_coarse, na.rm=TRUE))

message("trained ML model (GBM classifier) on ",
        length(unique(gbm_training_set$point_id)), " unique spatial points, ",
        nrow(gbm_training_set), " point-years")

message("GBM uses ", length(predictors), " MODIS predictors")

message("fit HMM model on ",
        length(unique(points_test$point_id)), " unique spatial points, ",
        nrow(points_test), " point-years")

message("HMM estimates of P:")
print(hmm_params_hat$P)  # Compare to pr_transition, pr_transition_predictions

message("Ground truth transition probabilities:")
print(pr_transition)

message("Transition probabilities in GBM predictions:")
print(pr_transition_predictions)

## Note: hmm_params_hat$pr_y and pr_Y_given_S are the transposes of the Upsilon matrix in the paper
message("HMM estimates of Pr[Y | S]:")
print(hmm_params_hat$pr_y)  # Compare to pr_Y_given_S

message("Pr[Y | S] computed on test set:")
print(pr_Y_given_S)
