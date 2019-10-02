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
                 make_option("--panel_length", default="long"))
opt <- parse_args(OptionParser(option_list=opt_list))
message("command line options: ", paste(sprintf("%s=%s", names(opt), opt), collapse=", "))
stopifnot(opt$landuse_set %in% c("binary", "ternary"))
stopifnot(opt$panel_length %in% c("long", "short"))

# TODO Should we experiment with a smaller GBM training set and larger HMM validation set?
training_fraction <- 0.50
stopifnot(0 < training_fraction && training_fraction < 1)

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
message("keeping ", nrow(points), " point-years, ", length(unique(points$point_id)), " unique spatial points")

point_id_train <- sample(unique(points$point_id), size=floor(training_fraction * length(unique(points$point_id))))

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

rf_test_confusion <- list()
gbm_test_confusion <- list()
pr_transition <- list()
pr_transition_predictions <- list()
hmm_params_hat <- list()

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

points_train <- subset(points, point_id %in% point_id_train)
points_test <- subset(points, !point_id %in% point_id_train)
table(points_train$validation_landuse_coarse)  # Careful, pasture is rare -- should make sure that training set is not 100% crops, same for test
table(points_test$validation_landuse_coarse)  # Similar comment about non-soy crops in ternary case

predictors <- c(evi_vars, nir_vars, mir_vars, red_vars, blue_vars)
model_formula <- reformulate(termlabels=predictors, response="validation_landuse_coarse")

model_rf <- randomForest(formula=model_formula, ntree=1000, nodesize=1, do.trace=TRUE,
                         mtry=floor(sqrt(length(predictors))), na.action=na.omit, data=points_train)
model_rf$confusion  # Large OOB error rate for pasture -- non-soy crops even worse in ternary case

points_test$landuse_predicted_rf <- predict(model_rf, type="response", newdata=points_test)
mean(is.na(points_test$landuse_predicted_rf))  # Around 0.11 -- missing when MODIS variables are NA
table(points_test$landuse_predicted_rf)
table(points_test$validation_landuse_coarse, points_test$landuse_predicted_rf)  # RF test confusion matrix, compare to OOB

## Fit GBM, see whether it beats RF -- takes around 10 minutes to fit
## GBM appears to beat RF for rare landuse in both binary and ternary case
gbm_filename <- sprintf("~/Dropbox/amazon_hmm_shared/data/gbm_model_for_embrapa_validation_%s_landuse_set_%s.rds",
                        opt$landuse_set, digest(points_train, algo="crc32"))
if(file.exists(gbm_filename)) {
    model_gbm <- readRDS(gbm_filename)
} else {
    model_gbm <- gbm(formula=model_formula, distribution="multinomial",  # Does this have to be bernoulli for two-class case?
                     data=subset(points, !is.na(validation_landuse_coarse)), n.trees=6500, interaction.depth=4, shrinkage=0.001, cv.folds=3)
    saveRDS(model_gbm, file=gbm_filename)  # Takes a while to fit -- save and reuse on future runs
}
ntrees_star <- gbm.perf(model_gbm, method="cv")  # Plot shows training and CV deviance as function of number of trees -- selects around 5000-6000 trees
message("GBM uses ", ntrees_star, " trees (tuning parameter selected by cross-validation)")

gbm_predictions <- predict(model_gbm, type="response", newdata=points_test, ntrees=ntrees_star)  # Outputs class probabilities
points_test$landuse_predicted_gbm <- colnames(gbm_predictions)[as.integer(apply(gbm_predictions, 1, which.max))]
mean(is.na(points_test$landuse_predicted_gbm))  # Zero (i.e. no missing predictions), unlike RF
table(points_test$landuse_predicted_gbm)
table(points_test$validation_landuse_coarse, points_test$landuse_predicted_gbm)  # GBM test confusion matrix, compare to RF
table(points_test$landuse_predicted_gbm, points_test$landuse_predicted_rf)  # GBM versus RF predictions

rf_test_confusion[[opt$landuse_set]] <- table(points_test$validation_landuse_coarse, points_test$landuse_predicted_rf)

## Careful, not diagonally dominant in ternary case (see non-soy crops entry), expect HMM to perform poorly
gbm_test_confusion[[opt$landuse_set]] <- table(points_test$validation_landuse_coarse, points_test$landuse_predicted_gbm)

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

pr_transition[[opt$landuse_set]] <- with(dtable, prop.table(table(validation_landuse_coarse, validation_landuse_coarse_next), 1))
pr_transition_predictions[[opt$landuse_set]] <- with(dtable, prop.table(table(landuse_predicted, landuse_predicted_next), 1))
message("pr_transition: ", paste(round(pr_transition[[opt$landuse_set]], 3), collapse=" "))
message("pr_transition_predictions: ", paste(round(pr_transition_predictions[[opt$landuse_set]], 3), collapse=" "))

list_of_initial_hmm_params <- get_initial_hmm_params(opt$landuse_set)

panel <- get_hmm_panel_from_points(dtable, discrete_y_varname="landuse_predicted")  # List of panel elements
baum_welch_time_homogeneous(panel[[1]], params=list_of_initial_hmm_params[[1]])  # Test, compare pi to panel[[1]]$y and panel[[1]]$validation_landuse

list_of_hmm_params_hat <- lapply(list_of_initial_hmm_params, function(initial_hmm_params) {
    return(em_parameter_estimates_time_homogeneous(panel=panel, params=initial_hmm_params, max_iter=50, epsilon=0.001))
})

hmm_params_hat[[opt$landuse_set]] <- list_of_hmm_params_hat[[which.max(sapply(list_of_hmm_params_hat, function(x) max(x$loglik)))]]  # Choose estimate with highest loglik
hmm_params_hat[[opt$landuse_set]]$P  # Compare to pr_transition, pr_transition_predictions
hmm_params_hat[[opt$landuse_set]]$pr_y  # Compare to with(points_test, prop.table(table(validation_landuse_coarse, landuse_predicted), margin=1))
message("hmm_params_hat$P: ", paste(round(hmm_params_hat[[opt$landuse_set]]$P, 3), collapse=" "))

## Run Viterbi and see whether test performance beats raw GBM
viterbi_paths <- lapply(panel, viterbi_path_time_homogeneous, params=hmm_params_hat[[opt$landuse_set]])
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
boots_filename <- sprintf("bootstrap_%s_landuse_set_%s_replications_%s.rds",
                          opt$landuse_set, n_boostrap_samples, digest(points_train, algo="crc32"))
if(file.exists(boots_filename)) {
    message("loading ", boots_filename)
    boots <- readRDS(boots_filename)
} else {
    message("did not find ", boots_filename, "; running now; time is ", Sys.time())
    boots <- replicate(n_boostrap_samples, run_bootstrap(), simplify=FALSE)  # TODO Slow, check whether saved file exists; if so, don't re-run
    saveRDS(boots, file=boots_filename)
}

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
    outfile <- sprintf("validation_long_panel_%s_bootstrap_%s_transition_hmm_and_gbm_predictions.png",
                       opt$landuse_set, landuse_underscore, landuse_underscore)
    ggsave(outfile, p, width=10, height=8)
    ## GBM test errors and HMM pr_y
    p <- (ggplot(df_boots, aes_string(x=test_error, y=hmm_error)) +
          geom_abline(slope=1, intercept=0, linetype=2, color="grey", size=1.2) +
          xlab(sprintf("error rate for land use predictions in bootstrapped test data, given true land use is %s", landuse)) +
          ylab("HMM estimate of error rate") +
          geom_point())
    outfile <- sprintf("validation_long_panel_%s_bootstrap_%s_misclassification_probability.png",
                       opt$landuse_set, landuse_underscore)
    ggsave(outfile, p, width=10, height=8)
}
outfile <- sprintf("validation_bootstrap_long_panel_%s_landuse_set_%s_replications.csv", opt$landuse_set, n_boostrap_samples)
write.csv(df_boots, file=outfile, row.names=FALSE)
