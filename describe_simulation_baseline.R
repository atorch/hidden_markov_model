library(data.table)
library(ggplot2)
library(optparse)
library(stargazer)
library(stringr)

source("hmm_functions.R")

opt_list <- list(make_option("--simulation_date", default="2019-10-21", type="character",
                             help="The date on which you ran simulation_baseline.R (the simulation output filenames include a date suffix)."))
opt <- parse_args(OptionParser(option_list=opt_list))

message("command line options: ", paste(sprintf("%s=%s", names(opt), opt), collapse=", "))

get_data_table_summarizing_single_simulation <- function(sim) {

    true_deforestation_prob <- sapply(sim$params$P_list,
                                      get_deforestation_prob_from_P)
    true_reforestation_prob <- sapply(sim$params$P_list,
                                      get_reforestation_prob_from_P)
    
    estimated_deforestation_prob_em <- sapply(sim$estimates$em_params_hat_best_likelihood$P_list,
                                              get_deforestation_prob_from_P)

    estimated_deforestation_prob_md <- sapply(sim$estimates$min_dist_params_hat_best_objfn$P_list,
                                              get_deforestation_prob_from_P)

    estimated_reforestation_prob_em <- sapply(sim$estimates$em_params_hat_best_likelihood$P_list,
                                                  get_reforestation_prob_from_P)

    estimated_reforestation_prob_md <- sapply(sim$estimates$min_dist_params_hat_best_objfn$P_list,
                                              get_reforestation_prob_from_P)

    muNaive <- mean(sapply(sim$simulation, function(vec) vec$y[1])==1) ## Observed forest cover in initial period
    
    P_hat_naive <- lapply(sim$estimates$M_Y_joint_hat, get_transition_probs_from_M_S_joint)
    estimate_deforestation_prob_naive <- sapply(P_hat_naive, get_deforestation_prob_from_P)
    estimate_reforestation_prob_naive <- sapply(P_hat_naive, get_reforestation_prob_from_P)
    
    ## EM means expectation-maximization, MD means minimum distance
    return(data.table(x=sim$deforestation_rates,
                      fixed_effect=sim$fixed_effect,
                      time=seq_along(sim$deforestation_rates),
                      sim_id=sim$id,
                      true_deforestation_probability=true_deforestation_prob,
                      true_reforestation_probability=true_reforestation_prob,
                      estimated_deforestation_probability_naive=estimate_deforestation_prob_naive,
                      estimated_deforestation_probability_em=estimated_deforestation_prob_em,
                      estimated_deforestation_probability_md=estimated_deforestation_prob_md,
                      estimated_reforestation_probability_naive=estimate_reforestation_prob_naive,
                      estimated_reforestation_probability_em=estimated_reforestation_prob_em,
                      estimated_reforestation_probability_md=estimated_reforestation_prob_md,
                      true_misclassification_probability_1=sim$params$pr_y[1,2],
                      true_misclassification_probability_2=sim$params$pr_y[2,1],
                      true_mu1 = sim$params$mu[1],
                      estimated_misclassification_probability_1_md=sim$estimates$min_dist_params_hat_best_objfn$pr_y[1,2],
                      estimated_misclassification_probability_1_em=sim$estimates$em_params_hat_best_likelihood$pr_y[1,2],
                      estimated_misclassification_probability_2_md=sim$estimates$min_dist_params_hat_best_objfn$pr_y[2,1],
                      estimated_misclassification_probability_2_em=sim$estimates$em_params_hat_best_likelihood$pr_y[2,1],
                      estimated_mu1_md=sim$estimates$min_dist_params_hat_best_objfn$mu[1],
                      estimated_mu1_naive=muNaive,
                      estimated_mu1_em=sim$estimates$em_params_hat_best_likelihood$mu[1])
           )
}

filename <- sprintf("output/baseline_simulation_%s_Desc.csv", opt$simulation_date)
message("Loading ", filename)
simulation_params <- fread(filename)
simulation_params$iter <- seq_len(nrow(simulation_params))

sim_dts <- lapply(seq_len(nrow(simulation_params)), function(iter) {
    filename <- sprintf("output/baseline_simulation_%s_iter_%s.rds", opt$simulation_date, iter)
    message("Loading ", filename)
    sim <- readRDS(filename)
    dt <- rbindlist(lapply(sim, get_data_table_summarizing_single_simulation))
    dt$iter <- iter
    return(dt)
})

dt <- rbindlist(sim_dts)

graph_colors  <- setNames(c("green","blue","purple"), c("Freq","MD","EM"))

id_vars <- c("time",
             "true_deforestation_probability",
             "true_reforestation_probability",
             "true_misclassification_probability_1",
             "true_misclassification_probability_2",
             "true_mu1",
             "x",
             "sim_id",
             "iter")
dt_melt  <- melt(dt, id.vars=id_vars)

dt_melt  <- merge(dt_melt, simulation_params, on="iter")

## This gives us the MSE, bias and variance for estimated_deforestation_probability_{naive,em,md}
dt_melt[prY11 == 90 &
        n_time_periods == 4 &
        time == 3 &
        defRtLast == 20 &
        n_points == 1000 &
        variable %like% 'deforestation_probability',
        list(mse = mean((value - true_deforestation_probability)^2),
             bias = mean(value - true_deforestation_probability),
             variance = var(value),
             N = .N),
        by = list(true_deforestation_probability, variable)]

n_points_disp_levels <- paste0("'N=",
                               str_trim(format(sort(unique(dt_melt$n_points)), big.mark=",",scientific=FALSE)),
                               " Points'")
dt_melt[, n_points_disp := factor(paste0("'N=",
                                         str_trim(format(n_points, big.mark=",",scientific=FALSE)),
                                         " Points'"),
                                  levels=n_points_disp_levels)]

dt_melt[variable %like% '_em', estimTypDisp := 'EM']
dt_melt[variable %like% '_md', estimTypDisp := 'MD']
dt_melt[variable %like% '_naive', estimTypDisp := 'Freq']

## Note: the order of the factor levels controls the order along the x-axis in several of the graphs
dt_melt[, estimTypDisp := factor(estimTypDisp, levels=c("Freq", "EM", "MD"))]

## TODO Could probably condition on iter instead of these multiple variables, but maybe this is more clear?
plt <- (ggplot(dt_melt[n_time_periods == 4 & prY11 == 90 & defRtLast == 20 & variable %like% 'deforestation_probability'],
               aes(y=value, x=estimTypDisp, group=variable)) +
        geom_boxplot() +
        geom_hline(aes(yintercept=true_deforestation_probability),linetype = 'dashed') +
        scale_y_continuous('Estimated Transition Probability',limits = c(0,.5)) +
        scale_x_discrete('Estimator') +
        facet_grid(paste0('P[',time,']') ~ n_points_disp,labeller = label_parsed) +
        theme_bw())
ggsave('output/deforestation_probability_different_sample_size.png', width=6, height=3, units='in')

plt <- (ggplot(dt_melt[n_points == 1000 & prY11 == 90 & defRtLast == 20 & variable %like% 'deforestation_probability'],
               aes(y=value, x = estimTypDisp, group=variable)) +
        geom_boxplot() +
        geom_hline(aes(yintercept=true_deforestation_probability),linetype = 'dashed')+
        scale_y_continuous('Estimated Transition Probability',limits = c(0,.5))+
        scale_x_discrete('Estimator') +
        facet_grid(paste0('P[',time,']')~paste0("'T=",n_time_periods," Periods'"),labeller = label_parsed)+
        theme_bw())
ggsave('output/deforestation_probability_different_nPeriods.png',width = 6,height = 4,units='in')

plt <- (ggplot(dt_melt[prY11 == 90 & n_time_periods == 4 & defRtLast == 20 & time == 1 & variable %like% 'misclassification_probability_1'],
               aes(y=value, x = estimTypDisp, group=variable)) +
        geom_boxplot() +
        geom_hline(aes(yintercept=true_misclassification_probability_1), linetype = 'dashed')+
        scale_x_discrete('Estimator', labels = function(x) str_replace_all(x,'estimated_misclassification_probability_1_','')) +
        scale_y_continuous('Estimated Misclassification Probability') +
        facet_wrap(~ n_points_disp,labeller = label_parsed) +
        theme_bw())
ggsave('output/misclassification_probability_different_sample_size.png', width = 6, height = 3, units='in')

plt <- (ggplot(dt_melt[prY11 == 90 & n_points == 1000 & defRtLast == 20 & time == 1 & variable %like% 'misclassification_probability_1'],
               aes(y=value, x = estimTypDisp, group=variable)) +
        geom_boxplot() +
        geom_hline(aes(yintercept=true_misclassification_probability_1), linetype = 'dashed')+
        scale_x_discrete('Estimator', labels = function(x) str_replace_all(x,'estimated_misclassification_probability_1_','')) +
        scale_y_continuous('Estimated Misclassification Probability')+
        facet_wrap(~paste0("T=",n_time_periods,' Periods'))+
        theme_bw())
ggsave('output/misclassification_probability_different_nPeriods.png',width = 6,height = 3,units='in')

plt <- (ggplot(dt_melt[n_time_periods == 4 & n_points == 1000 & defRtLast == 20 & time == 1 & variable %like% 'misclassification_probability_1'],
               aes(y=value, x = estimTypDisp, group=variable)) +
    geom_boxplot() +
    geom_hline(aes(yintercept=true_misclassification_probability_1), linetype = 'dashed')+
    scale_x_discrete('Estimator', labels = function(x) str_replace_all(x,'estimated_misclassification_probability_1_','')) +
    scale_y_continuous('Estimated Misclassification Probability')+
    facet_wrap(~paste0("Pr[Y=2|S=1]=", (1-as.integer(prY11) / 100)))+
    theme_bw())
ggsave('output/misclassification_probability_different_PrY11.png',width = 6,height = 4,units='in')


plt <- (ggplot(dt_melt[prY11 == 90 & time == 3 & n_time_periods == 4 & n_points == 1000 & variable %like% 'misclassification_probability_1'],
               aes(x=true_deforestation_probability, y=value, color = estimTypDisp,fill=estimTypDisp)) +
        geom_smooth() +
        geom_hline(aes(yintercept=true_misclassification_probability_1),linetype = 'dashed')+
        scale_color_manual('Estimator',values = graph_colors, labels = function(x) str_replace_all(x,'estimated_misclassification_probability_1_','')) +
        scale_fill_manual('Estimator',values = graph_colors, labels = function(x) str_replace_all(x,'estimated_misclassification_probability_1_','')) +
        scale_x_continuous(bquote("True Transition Probability(Pr["~ S[t==4]==2~"|"~S[t==3]==1~"])"))+
        scale_y_continuous('Estimated Misclassification Probability',limits = c(0,.3))+
        theme_bw())
ggsave('output/misclassification_probability_different_deforestation_rate.png', width = 6, height = 4,units = 'in')


plt <- (ggplot(dt_melt[ prY11 == 90 & n_points == 1000 & n_time_periods == 4 & variable %like% 'deforestation'],
               aes(x=defRtLast/100, y=value, color = estimTypDisp,fill=estimTypDisp)) +
        geom_smooth() +
        geom_line(aes(y=true_deforestation_probability),linetype='dashed',color='black')+
        scale_color_manual('Estimator',values = graph_colors, labels = function(x) str_replace_all(x,'estimated_deforestation_probability_','')) +
        scale_fill_manual('Estimator',values = graph_colors, labels = function(x) str_replace_all(x,'estimated_deforestation_probability_','')) +
        scale_x_continuous(bquote("True Transition Probability(Pr["~ S[t==4]==2~"|"~S[t==3]==1~"])"))+
        scale_y_continuous(bquote("Estimated Transition Probability"),limits = c(0,.4))+
        facet_grid(paste0('P[list(',time,')]')~.,label = label_parsed)+
        theme_bw())
ggsave('output/deforestation_probability_different_deforestation_rate.png', width = 6, height = 4,units = 'in')

plt <- (ggplot(dt_melt[ defRtLast==20 & time == 3 & n_time_periods == 4 & n_points == 1000 & variable %like% 'misclassification_probability_1'],
               aes(x=true_misclassification_probability_1, y=value, color = estimTypDisp,fill=estimTypDisp)) +
        geom_smooth() +
        geom_line(aes(y=true_misclassification_probability_1),linetype = 'dashed',color='black')+
        scale_color_manual('Estimator',values= graph_colors, labels = function(x) str_replace_all(x,'estimated_misclassification_probability_1_','')) +
        scale_fill_manual('Estimator',values = graph_colors, labels = function(x) str_replace_all(x,'estimated_misclassification_probability_1_','')) +
        scale_x_continuous("True Misclassification Probability(Pr[Y=2|S=1])")+
        scale_y_continuous('Estimated Misclassification Probability',limits = c(0,.4))+
        theme_bw())
ggsave('output/misclassification_probability_different_misclassification_rate.png', width = 6, height = 4,units = 'in')

plt <- (ggplot(dt_melt[defRtLast==20 & n_points == 1000 & n_time_periods == 4 & variable %like% 'deforestation'],
               aes(x=true_misclassification_probability_1, y=value, color = estimTypDisp,fill=estimTypDisp)) +
        geom_smooth() +
        geom_line(aes(y=true_deforestation_probability),linetype='dashed',color='black')+
        scale_color_manual('Estimator',values = graph_colors, labels = function(x) str_replace_all(x,'estimated_deforestation_probability_','')) +
        scale_fill_manual('Estimator', values = graph_colors, labels = function(x) str_replace_all(x,'estimated_deforestation_probability_','')) +
        scale_x_continuous("True Misclassification Probability(Pr[Y=2|S=1])")+
        scale_y_continuous("Estimated Transition Probability", limits = c(0,.4))+
        facet_grid(paste0('P[list(',time,')]')~., label = label_parsed)+
        theme_bw())
ggsave('output/deforestation_probability_different_misclassification_rate.png', width = 6, height = 4, units = 'in')
