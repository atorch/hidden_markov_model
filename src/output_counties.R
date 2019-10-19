library(data.table)
library(stargazer)
library(ggplot2)
library(viridis)
library(stringr)

source("hmm_functions.R")


fileReadVec <-c('2019-10-17 18:47:20')

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

county_df  <- NULL
for (f in seq_len(length(fileReadVec))){
    county_dfs <- list()
    for (i in seq_len(nrow(iterDat[fileID == f]))){
        counties  <- readRDS(paste0('county_simulation_',fileReadVec[f],'_iter_',i,'.rds'))
        county_dfs[[i]] <- rbindlist(lapply(counties, get_data_table_summarizing_single_county_simulation))
    }
    county_dff <- rbindlist(county_dfs,idcol = 'simID')[,fileID := f]
    if (exists('county_df')){
        county_df  <- rbind(county_df,county_dff,fill=TRUE,use.names=TRUE)
    }else{
        county_df  <- county_dff
    }
}

##colors for graph
graphCol  <- setNames(c('green','blue','purple'),c('Frequency','MD','EM'))

##Make Graphs
county_df_melt  <- melt(county_df,id.vars = c('time','county_id',
                                               'true_deforestation_probability','true_reforestation_probability','true_misclassification_probability_1','true_misclassification_probability_2','true_mu1','x','fixed_effect','simID','fileID'))

county_df_melt  <- merge(county_df_melt, iterDat)

##Elimate Duplicated Rows
county_df_melt  <- county_df_melt[!iterDat[duplicated(iterDat[,2:10]),list(simID,fileID)],on=c('simID','fileID')]

county_df_melt[ prY11 == 90 & n_time_periods == 4 & time == 3 & defRtLast == 20 & n_points_per_county == 1000 & variable %like% 'deforestation_probability', list(mse = mean((value - true_deforestation_probability)^2),
                                                                  bias = mean(value - true_deforestation_probability),
                                                                  variance = var(value),
                                                                  N = .N),
               by = list(true_deforestation_probability,variable)]

county_df_melt[,n_points_per_county_disp := factor(paste0("'N=",n_points_per_county, " Points'"), levels = paste0("'N=",sort(unique(n_points_per_county)), " Points'"))]

county_df_melt[variable %like% '_em', estimTypDisp := 'EM']
county_df_melt[variable %like% '_md', estimTypDisp := 'MD']
county_df_melt[variable %like% '_naive', estimTypDisp := 'Frequency']
##Graphs
plt <- ggplot(county_df_melt[ n_time_periods == 4 & prY11 == 90 & defRtLast == 20 & variable %like% 'deforestation_probability'], aes( y=value, x = estimTypDisp, group=variable)) +
    geom_boxplot() +
    geom_hline(aes(yintercept=true_deforestation_probability),linetype = 'dashed')+
    scale_y_continuous('Estimated Transition Probability',limits = c(0,.5))+
    scale_x_discrete('Estimator') +
    facet_grid(paste0('P[',time,']')~n_points_per_county_disp,labeller = label_parsed)+
    theme_bw()
ggsave('deforestation_probability_different_sample_size.png',width = 6,height=3,units='in')

plt <- ggplot(county_df_melt[ n_points_per_county == 1000 & prY11 == 90 & defRtLast == 20 & variable %like% 'deforestation_probability'], aes( y=value, x = estimTypDisp, group=variable)) +
    geom_boxplot() +
    geom_hline(aes(yintercept=true_deforestation_probability),linetype = 'dashed')+
    scale_y_continuous('Estimated Transition Probability',limits = c(0,.5))+
    scale_x_discrete('Estimator') +
    facet_grid(paste0('P[',time,']')~paste0("'T=",n_time_periods," Periods'"),labeller = label_parsed)+
    theme_bw()
ggsave('deforestation_probability_different_nPeriods.png',width = 6,height = 4,units='in')


plt <- ggplot(county_df_melt[prY11 == 90 &  n_time_periods == 4 & defRtLast == 20 & time == 1 & variable %like% 'misclassification_probability_1'], aes( y=value, x = estimTypDisp, group=variable)) +
    geom_boxplot() +
    geom_hline(aes(yintercept=true_misclassification_probability_1), linetype = 'dashed')+
    scale_x_discrete('Estimator', labels = function(x) str_replace_all(x,'estimated_misclassification_probability_1_','')) +
    scale_y_continuous('Estimated Misclassification Probability')+
    facet_wrap(~n_points_per_county_disp,labeller = label_parsed)+
    theme_bw()
ggsave('misclassification_probability_different_sample_size.png',width = 6,height = 3,units='in')

plt <- ggplot(county_df_melt[prY11 == 90 &  n_points_per_county == 1000 & defRtLast == 20 & time == 1 & variable %like% 'misclassification_probability_1'], aes( y=value, x = estimTypDisp, group=variable)) +
    geom_boxplot() +
    geom_hline(aes(yintercept=true_misclassification_probability_1), linetype = 'dashed')+
    scale_x_discrete('Estimator', labels = function(x) str_replace_all(x,'estimated_misclassification_probability_1_','')) +
    scale_y_continuous('Estimated Misclassification Probability')+
    facet_wrap(~paste0("T=",n_time_periods,' Periods'))+
    theme_bw()
ggsave('misclassification_probability_different_nPeriods.png',width = 6,height = 3,units='in')

plt <- ggplot(county_df_melt[n_time_periods == 4 &  n_points_per_county == 1000 & defRtLast == 20 & time == 1 & variable %like% 'misclassification_probability_1'], aes( y=value, x = estimTypDisp, group=variable)) +
    geom_boxplot() +
    geom_hline(aes(yintercept=true_misclassification_probability_1), linetype = 'dashed')+
    scale_x_discrete('Estimator', labels = function(x) str_replace_all(x,'estimated_misclassification_probability_1_','')) +
    scale_y_continuous('Estimated Misclassification Probability')+
    facet_wrap(~paste0("Pr[Y=2|S=1]=",(1-as.integer(prY11)/100)))+
    theme_bw()
ggsave('misclassification_probability_different_PrY11.png',width = 6,height = 4,units='in')


plt <- ggplot(county_df_melt[ prY11 == 90 & time == 3 & n_time_periods == 4 & n_points_per_county == 1000 & variable %like% 'misclassification_probability_1'], aes(x=true_deforestation_probability, y=value, color = estimTypDisp,fill=estimTypDisp)) +
    geom_smooth() +
    geom_hline(aes(yintercept=true_misclassification_probability_1),linetype = 'dashed')+
    scale_color_manual('Estimator',values = graphCol, labels = function(x) str_replace_all(x,'estimated_misclassification_probability_1_','')) +
    scale_fill_manual('Estimator',values = graphCol, labels = function(x) str_replace_all(x,'estimated_misclassification_probability_1_','')) +
    scale_x_continuous(bquote("True Transition Probability(Pr["~ S[t==4]==2~"|"~S[t==3]==1~"])"))+
    scale_y_continuous('Estimated Misclassification Probability',limits = c(0,.3))+
    theme_bw()
ggsave('misclassification_probability_different_deforestation_rate.png', width = 6, height = 4,units = 'in')


plt <- ggplot(county_df_melt[ prY11 == 90 &  n_points_per_county == 1000 & n_time_periods == 4 & variable %like% 'deforestation'], aes(x=defRtLast/100, y=value, color = estimTypDisp,fill=estimTypDisp)) +
    geom_smooth() +
    geom_line(aes(y=true_deforestation_probability),linetype='dashed',color='black')+
    scale_color_manual('Estimator',values = graphCol, labels = function(x) str_replace_all(x,'estimated_deforestation_probability_','')) +
    scale_fill_manual('Estimator',values = graphCol, labels = function(x) str_replace_all(x,'estimated_deforestation_probability_','')) +
    scale_x_continuous(bquote("True Transition Probability(Pr["~ S[t==4]==2~"|"~S[t==3]==1~"])"))+
    scale_y_continuous(bquote("Estimated Transition Probability"),limits = c(0,.4))+
    facet_grid(paste0('P[list(',time,')]')~.,label = label_parsed)+
    theme_bw()
ggsave('deforestation_probability_different_deforestation_rate.png', width = 6, height = 4,units = 'in')

plt <- ggplot(county_df_melt[ defRtLast==20 & time == 3 & n_time_periods == 4 & n_points_per_county == 1000 & variable %like% 'misclassification_probability_1'], aes(x=true_misclassification_probability_1, y=value, color = estimTypDisp,fill=estimTypDisp)) +
    geom_smooth() +
    geom_line(aes(y=true_misclassification_probability_1),linetype = 'dashed',color='black')+
    scale_color_manual('Estimator',values= graphCol, labels = function(x) str_replace_all(x,'estimated_misclassification_probability_1_','')) +
    scale_fill_manual('Estimator',values = graphCol, labels = function(x) str_replace_all(x,'estimated_misclassification_probability_1_','')) +
    scale_x_continuous("True Misclassification Probability(Pr[Y=2|S=1])")+
    scale_y_continuous('Estimated Misclassification Probability',limits = c(0,.4))+
    theme_bw()
ggsave('misclassification_probability_different_misclassification_rate.png', width = 6, height = 4,units = 'in')


plt <- ggplot(county_df_melt[ defRtLast==20& n_points_per_county == 1000 & n_time_periods == 4 & variable %like% 'deforestation'], aes(x=true_misclassification_probability_1, y=value, color = estimTypDisp,fill=estimTypDisp)) +
    geom_smooth() +
    geom_line(aes(y=true_deforestation_probability),linetype='dashed',color='black')+
    scale_color_manual('Estimator',values = graphCol, labels = function(x) str_replace_all(x,'estimated_deforestation_probability_','')) +
    scale_fill_manual('Estimator', values = graphCol, labels = function(x) str_replace_all(x,'estimated_deforestation_probability_','')) +
    scale_x_continuous("True Misclassification Probability(Pr[Y=2|S=1])")+
    scale_y_continuous("Estimated Transition Probability",limits = c(0,.4))+
    facet_grid(paste0('P[list(',time,')]')~.,label = label_parsed)+
    theme_bw()
ggsave('deforestation_probability_different_misclassification_rate.png', width = 6, height = 4, units = 'in')


mseFunc <- function(trueVal,keyword){
    prDat <- county_df_melt[ variable %like% keyword, list(rmse = sqrt(mean((value - get(trueVal))^2)),
                                                                                  bias = mean(value - get(trueVal)),
                                                                                  sd = sd(value),
                                                                                  N = .N),
                               by = list(variable,fileID,simID,time)]
    return(prDat)
}


nVec <- c(100,500,1000)

colsT <- iterDat[n_points_per_county %in% nVec & n_counties == 100  & mu1==90 & defRt1==4&defRtMid==10 & defRtLast==20 & prY11==90 & prY22==80 & n_points_per_county == 1000]
colsMrgT <- colsT[,list(simID,fileID,nPts = n_points_per_county)]

defForestPrDatT <- merge(mseFunc('true_deforestation_probability','deforestation_probability'),colsMrgT)
miscPr1DatT <- merge(mseFunc('true_misclassification_probability_1','misclassification_probability_1'),colsMrgT)


cols <- iterDat[n_points_per_county %in% nVec & n_counties == 100 & n_time_periods==4 & mu1==90 & defRt1==4&defRtMid==10 & defRtLast==20 & prY11==90 & prY22==80]
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
nmVec <- c('P_{S_{1}}=.9','\\Upsilon(2,1)=.1','\\Upsilon(1,2)=.2','P_{1}(1,2)=.04','P_{1}(2,1)=.02','P_{2}(1,2)=.1','P_{2}(2,1)=.02','P_{3}(1,2)=.2','P_{3}(2,1)=.02')

sink('mcTable.tex')
cat('\\begin{tabular}{rr@{\\hskip .3in}ccc@{\\hskip .4in}ccc@{\\hskip .4in}ccc}\n')
cat('\\hline\n')
cat('& & ',paste0('\\multicolumn{3}{c}{N=',nVec,'}',collapse='&'),'\\\\\n')
cat('\\hline\n')
cat('&  ',rep('&Freq & MD & EM',3),'\\\\\n')
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
                       
