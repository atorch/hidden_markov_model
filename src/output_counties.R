library(data.table)
library(stargazer)
library(ggplot2)
library(viridis)
library(stringr)

source("hmm_functions.R")


fileReadVec <-c('2019-09-24 19:24:42', '2019-09-25 21:50:24')

get_data_table_summarizing_single_county_simulation <- function(county) {

    true_deforestation_prob <- sapply(county$params$P_list,
                                      get_deforestation_prob_from_P)
    true_reforestation_prob <- sapply(county$params$P_list,
                                      get_reforestation_prob_from_P)
    
    estimated_deforestation_prob_em <- sapply(county$estimates$em_params_hat_best_likelihood$P_list,
                                              get_deforestation_prob_from_P)

    estimated_deforestation_prob_md <- sapply(county$estimates$min_dist_params_hat_best_objfn$P_list,
                                              get_deforestation_prob_from_P)

    estimated_reforestation_prob_em <- sapply(county$estimates$em_params_hat_best_likelihood$P_list,
                                                  get_reforestation_prob_from_P)

    estimated_reforestation_prob_md <- sapply(county$estimates$min_dist_params_hat_best_objfn$P_list,
                                              get_reforestation_prob_from_P)

    muNaive <- mean(sapply(county$simulation, function(vec) vec$y[1])==1) ##Observed forest cover in initial period
    
    P_hat_naive <- lapply(county$estimates$M_Y_joint_hat, get_transition_probs_from_M_S_joint)
    estimate_deforestation_prob_naive <- sapply(P_hat_naive, get_deforestation_prob_from_P)
    estimate_reforestation_prob_naive <- sapply(P_hat_naive, get_reforestation_prob_from_P)

    
    ## EM means expectation-maximization, MD means minimum distance
    return(data.table(x=county$X,
                      fixed_effect=county$fixed_effect,
                      time=seq_along(county$X),
                      county_id=county$id,

                      
                      true_deforestation_probability=true_deforestation_prob,
                      true_reforestation_probability=true_reforestation_prob,
                      estimated_deforestation_probability_naive=estimate_deforestation_prob_naive,
                      estimated_deforestation_probability_em=estimated_deforestation_prob_em,
                      estimated_deforestation_probability_md=estimated_deforestation_prob_md,
                      estimated_reforestation_probability_naive=estimate_reforestation_prob_naive,
                      estimated_reforestation_probability_em=estimated_reforestation_prob_em,
                      estimated_reforestation_probability_md=estimated_reforestation_prob_md,

                      true_misclassification_probability_1  = county$params$pr_y[1,2],
                      true_misclassification_probability_2  = county$params$pr_y[2,1],
                      true_mu1 = county$params$mu[1],
                      estimated_misclassification_probability_1_md  = county$estimates$min_dist_params_hat_best_objfn$pr_y[1,2],
                      estimated_misclassification_probability_1_em  = county$estimates$em_params_hat_best_likelihood$pr_y[1,2],
                      estimated_misclassification_probability_2_md  = county$estimates$min_dist_params_hat_best_objfn$pr_y[2,1],
                      estimated_misclassification_probability_2_em  = county$estimates$em_params_hat_best_likelihood$pr_y[2,1],
                      estimated_mu1_md  = county$estimates$min_dist_params_hat_best_objfn$mu[1],
                      estimated_mu1_naive  = muNaive,
                      estimated_mu1_em  = county$estimates$em_params_hat_best_likelihood$mu[1]
                      )
           )
}


##Import data
iterDat <- rbindlist(lapply(fileReadVec, function(x) fread(paste0('county_simulation_',x,'_Desc.csv'))[,simID := .I]),idcol='fileID')
iterDat <- iterDat[!(simID == 7 & fileID == 1)] ##TEMPORARY FIX
county_df  <- NULL
for (f in seq_len(length(fileReadVec))){
    county_dfs <- list()
    for (i in seq_len(nrow(iterDat[fileID == f]))){
        counties  <- readRDS(paste0('county_simulation_',fileReadVec[f],'_iter_',i,'.rds'))
        county_dfs[[i]] <- rbindlist(lapply(counties, get_data_table_summarizing_single_county_simulation))
    }
    county_dff <- rbindlist(county_dfs,idcol = 'simID')[,fileID := f]
    if (exists('county_df')){
        county_df  <- rbind(county_df,county_dff)
    }else{
        county_df  <- county_dff
    }
}

##Make Graphs
county_df_melt  <- melt(county_df,id.vars = c('time','county_id',
                                               'true_deforestation_probability','true_reforestation_probability','true_misclassification_probability_1','true_misclassification_probability_2','true_mu1','x','fixed_effect','simID','fileID'))

county_df_melt  <- merge(county_df_melt, iterDat)

county_df_melt[ prY11 == 90 & n_time_periods == 4 & time == 3 & defRtLast == 20 & n_points_per_county == 1000 & variable %like% 'deforestation_probability', list(mse = mean((value - true_deforestation_probability)^2),
                                                                  bias = mean(value - true_deforestation_probability),
                                                                  variance = var(value),
                                                                  N = .N),
               by = list(true_deforestation_probability,variable)]

county_df_melt[,n_points_per_county_disp := factor(paste0("'",n_points_per_county, " Points'"), levels = paste0("'",sort(unique(n_points_per_county)), " Points'"))]

county_df_melt[variable %like% '_em', estimTypDisp := 'EM']
county_df_melt[variable %like% '_md', estimTypDisp := 'MD']
county_df_melt[variable %like% '_naive', estimTypDisp := 'Naive']
##Graphs
plt <- ggplot(county_df_melt[ n_time_periods == 4 & prY11 == 90 & defRtLast == 20 & variable %like% 'deforestation_probability'], aes( y=value, x = estimTypDisp, group=variable)) +
    geom_boxplot() +
    geom_hline(aes(yintercept=true_deforestation_probability),linetype = 'dashed')+
    ggtitle('Deforestation Probabilities') +
    scale_y_continuous('Estimated Deforestation Probability',limits = c(0,.5))+
    scale_x_discrete('Estimator') +
    facet_grid(paste0('P[list(',time,',',time+1,')]')~n_points_per_county_disp,labeller = label_parsed)+
    theme_bw()
ggsave('deforestation_probability_different_sample_size.png',width = 8,height=5,units='in')

plt <- ggplot(county_df_melt[ n_points_per_county == 1000 & prY11 == 90 & defRtLast == 20 & variable %like% 'deforestation_probability'], aes( y=value, x = estimTypDisp, group=variable)) +
    geom_boxplot() +
    geom_hline(aes(yintercept=true_deforestation_probability),linetype = 'dashed')+
    ggtitle('Deforestation Probabilities') +
    scale_y_continuous('Estimated Deforestation Probability',limits = c(0,.5))+
    scale_x_discrete('Estimator') +
    facet_grid(paste0('P[list(',time,',',time+1,')]')~paste0("'N Periods='~",n_time_periods),labeller = label_parsed)+
    theme_bw()
ggsave('deforestation_probability_different_nPeriods.png',width = 8,height=5,units='in')


plt <- ggplot(county_df_melt[prY11 == 90 &  n_time_periods == 4 & defRtLast == 20 & time == 1 & variable %like% 'misclassification_probability_1'], aes( y=value, x = estimTypDisp, group=variable)) +
    geom_boxplot() +
    geom_hline(aes(yintercept=true_misclassification_probability_1), linetype = 'dashed')+
    ggtitle('Misclassification Probabilities') +
    scale_x_discrete('Estimator', labels = function(x) str_replace_all(x,'estimated_misclassification_probability_1_','')) +
    scale_y_continuous('Estimated Misclassification Probability')+
    facet_wrap(~n_points_per_county_disp,labeller = label_parsed)+
    theme_bw()
ggsave('misclassification_probability_different_sample_size.png',width = 8,height=5,units='in')

plt <- ggplot(county_df_melt[prY11 == 90 &  n_points_per_county == 1000 & defRtLast == 20 & time == 1 & variable %like% 'misclassification_probability_1'], aes( y=value, x = estimTypDisp, group=variable)) +
    geom_boxplot() +
    geom_hline(aes(yintercept=true_misclassification_probability_1), linetype = 'dashed')+
    ggtitle('Misclassification Probabilities') +
    scale_x_discrete('Estimator', labels = function(x) str_replace_all(x,'estimated_misclassification_probability_1_','')) +
    scale_y_continuous('Estimated Misclassification Probability')+
    facet_wrap(~paste0(n_time_periods,' Periods'))+
    theme_bw()
ggsave('misclassification_probability_different_nPeriods.png',width = 8,height=5,units='in')

plt <- ggplot(county_df_melt[n_time_periods == 4 &  n_points_per_county == 1000 & defRtLast == 20 & time == 1 & variable %like% 'misclassification_probability_1'], aes( y=value, x = estimTypDisp, group=variable)) +
    geom_boxplot() +
    geom_hline(aes(yintercept=true_misclassification_probability_1), linetype = 'dashed')+
    ggtitle('Misclassification Probabilities') +
    scale_x_discrete('Estimator', labels = function(x) str_replace_all(x,'estimated_misclassification_probability_1_','')) +
    scale_y_continuous('Estimated Misclassification Probability')+
    facet_wrap(~paste0('Y[list(1,1)]==',(1-as.integer(prY11)/100)),labeller = label_parsed)+
    theme_bw()
ggsave('misclassification_probability_different_PrY11.png',width = 8,height=5,units='in')


plt <- ggplot(county_df_melt[ prY11 == 90 & time == 3 & n_time_periods == 4 & n_points_per_county == 1000 & variable %like% 'misclassification_probability_1'], aes(x=true_deforestation_probability, y=value, color = estimTypDisp,fill=estimTypDisp)) +
    geom_smooth() +
    geom_hline(aes(yintercept=true_misclassification_probability_1),linetype = 'dashed')+
    ggtitle('Misclassification Probabilities as Deforestation Rates Change') +
    scale_color_viridis('Estimator',discrete = TRUE, labels = function(x) str_replace_all(x,'estimated_misclassification_probability_1_','')) +
    scale_fill_viridis('Estimator',discrete = TRUE, labels = function(x) str_replace_all(x,'estimated_misclassification_probability_1_','')) +
    scale_x_continuous(bquote("True Deforestation Probability("~ P[list(3,4)]~")"))+
    scale_y_continuous('Estimated Misclassification Probability')+
    theme_bw()
ggsave('misclassification_probability_different_deforestation_rate.png', width = 8, height = 5)


plt <- ggplot(county_df_melt[ prY11 == 90 &  n_points_per_county == 1000 & n_time_periods == 4 & variable %like% 'deforestation'], aes(x=defRtLast/100, y=value, color = estimTypDisp,fill=estimTypDisp)) +
    geom_smooth() +
    geom_line(aes(y=true_deforestation_probability),linetype='dashed',color='black')+
    ggtitle('Deforestation Probability Estimator Accuracy') +
    scale_color_viridis('Estimator',discrete = TRUE, labels = function(x) str_replace_all(x,'estimated_deforestation_probability_','')) +
    scale_fill_viridis('Estimator',discrete = TRUE, labels = function(x) str_replace_all(x,'estimated_deforestation_probability_','')) +
    scale_x_continuous(bquote("True Deforestation Probability("~ P[list(3,4)]~")"))+
    scale_y_continuous(bquote("Estimated Deforestation Probability("~ P[list(3,4)]~")"))+
    facet_wrap(~paste0('P[list(',time,',',time+1,')]'),label = label_parsed)+
    theme_bw()
ggsave('deforestation_probability_different_deforestation_rate.png', width = 8, height = 5)

mseFunc <- function(trueVal,keyword){
    prDat <- county_df_melt[ variable %like% keyword, list(rmse = sqrt(mean((value - get(trueVal))^2)),
                                                                                  bias = mean(value - get(trueVal)),
                                                                                  sd = sd(value),
                                                                                  N = .N),
                               by = list(variable,fileID,simID,time)]
    return(prDat)
}


nVec <- c(100,500,1000)
cols <- iterDat[n_points_per_county %in% nVec & n_counties == 50 & n_time_periods==4 & mu1==90 & defRt1==4&defRtMid==10 & defRtLast==20 & prY11==90 & prY22==80]
colsMrg <- cols[,list(simID,fileID,nPts = n_points_per_county)]


defForestPrDat <- merge(mseFunc('true_deforestation_probability','deforestation_probability'),colsMrg)
reForestPrDat <- merge(mseFunc('true_reforestation_probability','reforestation_probability'),colsMrg)
miscPr1Dat <- merge(mseFunc('true_misclassification_probability_1','misclassification_probability_1'),colsMrg)
miscPr2Dat <- merge(mseFunc('true_misclassification_probability_2','misclassification_probability_2'), colsMrg)
muPrDat <- merge(mseFunc('true_mu1','mu1'), colsMrg)



outFunc <- function(dat,nPtVar,varTyp,errTyp,ti=1){
    a <- dat[variable %like% varTyp  & time==ti & nPts == nPtVar,sprintf('%.3f',get(errTyp))]
    if (length(a) ==0){
        return('')
    }else{
        return(paste0(a))
    }
}



##Make output table 
matVec <- c('muPrDat','miscPr1Dat','miscPr2Dat',rep(c('defForestPrDat','reForestPrDat'),3))
tVec <- c(1,1,1,1,1,2,2,3,3)
nmVec <- c('\\mu=.9','\\Upsilon(1,2)=.1','\\Upsilon(2,1)=.2','P_{1,2}(1,2)=.04','P_{1,2}(2,1)=.02','P_{2,3}(1,2)=.1','P_{2,3}(2,1)=.02','P_{3,4}(1,2)=.2','P_{3,4}(2,1)=.02')

sink('mcTable.tex')
cat('\\begin{tabular}{rr@{\\hskip .3in}ccc@{\\hskip .4in}ccc@{\\hskip .4in}ccc}\n')
cat('\\hline\n')
cat('& & ',paste0('\\multicolumn{3}{c}{N=',nVec,'}',collapse='&'),'\\\\\n')
cat('\\hline\n')
cat('&  ',rep('&Naive & MD & EM',3),'\\\\\n')
cat('\\hline\n')
for(j in 1:length(matVec)){
    cat('&Bias &', paste0(sapply(c(100,500,1000),function(n) paste0(sapply(c('naive','md','em'),
                                                                           function(e) outFunc(get(matVec[j]),n,e,'bias',tVec[j])),collapse='&')),collapse='&'),'\\\\\n')
    cat('$',nmVec[j],'$& s.d. &', paste0(sapply(c(100,500,1000),function(n) paste0(sapply(c('naive','md','em'),
                                                                                          function(e) outFunc(get(matVec[j]),n,e,'sd',tVec[j])),collapse='&')),collapse='&'),'\\\\\n')
    cat('&RMSE &', paste0(sapply(c(100,500,1000),function(n) paste0(sapply(c('naive','md','em'),
                                                                           function(e) outFunc(get(matVec[j]),n,e,'rmse',tVec[j])),collapse='&')),collapse='&'),'\\\\\n')
    cat('\\\\\\\\')
}
cat('\\end{tabular}')
sink()


## cat('&Bias &', paste0(sapply(c(100,500,1000),function(n) paste0(sapply(c('naive','md','em'),
##                                                                        function(e) outFunc(miscPr2Dat,n,e,'bias')),collapse='&')),collapse='&'),'\\\\\n')
## cat('$\\Upsilon(2,1)=.2$& s.d. &', paste0(sapply(c(100,500,1000),function(n) paste0(sapply(c('naive','md','em'),
##                                                                                            function(e) outFunc(miscPr2Dat,n,e,'sd')),collapse='&')),collapse='&'),'\\\\\n')
## cat('&RMSE &', paste0(sapply(c(100,500,1000),function(n) paste0(sapply(c('naive','md','em'),
##                                                                        function(e) outFunc(miscPr2Dat,n,e,'rmse')),collapse='&')),collapse='&'),'\\\\\n')
## cat('\\\\\\\\')
## cat('&Bias &', paste0(sapply(c(100,500,1000),function(n) paste0(sapply(c('naive','md','em'),
##                                                                        function(e) outFunc(defForestPrDat,n,e,'bias',1)),collapse='&')),collapse='&'),'\\\\\n')
## cat('$P_{1}(1,2)=.04$& s.d. &', paste0(sapply(c(100,500,1000),function(n) paste0(sapply(c('naive','md','em'),
##                                                                                            function(e) outFunc(defForestPrDat,n,e,'sd',1)),collapse='&')),collapse='&'),'\\\\\n')
## cat('&RMSE &', paste0(sapply(c(100,500,1000),function(n) paste0(sapply(c('naive','md','em'),
##                                                                        function(e) outFunc(defForestPrDat,n,e,'rmse',1)),collapse='&')),collapse='&'),'\\\\\n')
## cat('\\\\\\\\')
## cat('&Bias &', paste0(sapply(c(100,500,1000),function(n) paste0(sapply(c('naive','md','em'),
##                                                                        function(e) outFunc(defForestPrDat,n,e,'bias',1)),collapse='&')),collapse='&'),'\\\\\n')
## cat('$P_{1}(2,1)=.02$& s.d. &', paste0(sapply(c(100,500,1000),function(n) paste0(sapply(c('naive','md','em'),
##                                                                                            function(e) outFunc(defForestPrDat,n,e,'sd',1)),collapse='&')),collapse='&'),'\\\\\n')
## cat('&RMSE &', paste0(sapply(c(100,500,1000),function(n) paste0(sapply(c('naive','md','em'),
##                                                                        function(e) outFunc(defForestPrDat,n,e,'rmse',1)),collapse='&')),collapse='&'),'\\\\\n')


## stargazer(simulation_df,
##           summary=FALSE,
##           rownames=FALSE,
##           out=simulation_df_outfile)
                       
