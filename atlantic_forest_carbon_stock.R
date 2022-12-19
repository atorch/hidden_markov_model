## Preliminaries
library(ggplot2)
library(data.table)
library(stargazer)
library(scales)
library(stringr)
library(pbapply)
library(optparse)

rm(list=ls())

opt_list <- list(make_option("--input_file", default="estimated_deforestation_rates_2022_01_14_use_md_as_initial_values_for_em.csv"))

opt <- parse_args(OptionParser(option_list=opt_list))

##input file
mapbiomasDat <- fread(opt$input_file)

carbonToCO2 <- 44/12 # (44 units CO2/12 units C) from https://www.epa.gov/energy/greenhouse-gases-equivalencies-calculator-calculations-and-references
tonsCO2ToDollars2040 <- 103 #from https://www.whitehouse.gov/wp-content/uploads/2021/02/TechnicalSupportDocument_SocialCostofCarbonMethaneNitrousOxide.pdf?source=email Using 2.5% discount rate for 2040 emissions
tonsCO2ToDollars2020 <- 76 #from https://www.whitehouse.gov/wp-content/uploads/2021/02/TechnicalSupportDocument_SocialCostofCarbonMethaneNitrousOxide.pdf?source=email Using 2.5% discount rate for 2040 emissions
carbonToDollars2020 <- carbonToCO2 * tonsCO2ToDollars2020
carbonToDollars2040 <- carbonToCO2 * tonsCO2ToDollars2040  
hectaresAtlanticForest <- 142274200 #From https://mapbiomas.org/en/mapbiomas-launches-1
baseYear <- 2010


#Colorblind friendly palettes http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# The palette with black:
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

carbonStockResDir <- './carbon_stock_results'
carbonStockFiles <- list.files(carbonStockResDir,full.names = TRUE)

##Data process
restrictedDat <- mapbiomasDat[ pr_y_diag_dominant_ml == TRUE &  n_mapbiomas_classes  %in% c(2,3)]

tilesWithViterbi <- str_match(carbonStockFiles,'Carbon_([1-9][0-9]+1\\_[1-9][0-9]+1)')[,2]
setdiff(restrictedDat[,unique(window_id)],tilesWithViterbi)

## List tiles not in restricted dat file, but for which we have a viterbi file
tilesInCarbonFiles <- str_match(list.files('./carbon_stock_results'),'Carbon_([1-9][0-9]+_[1-9][0-9]+)')[,2]
setdiff(tilesInCarbonFiles,restrictedDat[,unique(window_id)])



##Temporary reforestation values
#restrictedDat[,':='(reforestation_rate_md = NULL, reforestation_rate_ml=NULL)]
setnames(restrictedDat,
         paste('pr_agriculture_and_pasture_to_forest',c('md','ml','freq'),sep='_'),
         paste('reforestation_rateAP',c('md','ml','freq'),sep='_'))
         
 
## ##Collapse data for trend graphs
meltTrendDat <- melt(restrictedDat, id.vars = c('year','window_id'),measure = patterns('deforestation_rate|reforestation_rate|fraction_forest'))
meltTrendDat[,estimator := str_split_fixed(variable,'_',3)[,3]]
meltTrendDat[,measurement := paste(str_split_fixed(variable,'_',3)[,1],
                                   str_split_fixed(variable,'_',3)[,2],sep='_')]
meltTrendDat[measurement == 'fraction_forest', fraction_forest_tmp := value]
meltTrendDat[, fraction_forest  := min(fraction_forest_tmp,na.rm=TRUE),by = list(year,window_id,estimator)]
meltTrendDat[measurement == 'deforestation_rate',weights := fraction_forest, by = list(year,window_id,estimator)]
meltTrendDat[measurement == 'reforestation_rate',weights := 1-min(fraction_forest,na.rm=TRUE), by = list(year,window_id,estimator)] 
meltTrendDat[measurement == 'fraction_forest',weights := 1, by = list(year,window_id,estimator)] 

##Rates of Change in forest using HMM
meltTrendDat[,measurement := factor(measurement,levels = c('deforestation_rate','reforestation_rate','fraction_forest'))]
meltTrendDatGraph <- meltTrendDat[,list(value=weighted.mean(value,weights)),
                                  by = list(year,estimator,measurement)]

ggplot(meltTrendDatGraph[measurement!='reforestation_rateAP' & estimator !='md'],aes(y=value, x= year, color = estimator ))+geom_point()+
    stat_smooth() +
    theme_bw() +
    scale_y_continuous('Rate', labels = percent) +
    scale_color_manual('Estimator',values = cbPalette,breaks = c('ml','freq'),labels = c('HMM','Raw'))+
    facet_grid(measurement~., scales='free_y',
               labeller = as_labeller(c('deforestation_rate' = 'Deforestation Rate','fraction_forest' = 'Fraction Forest','reforestation_rate' = 'Reforestation Rate'))) +
    labs(caption = 'Points reflect the HMM parameters aggregated over all of the tiles\n(where the tile-specific deforestation rates are weighted by the fraction forest and the reforestation\nrates are weighted by the fraction not-forest.) Lines reflect a Loess trend.')+
    theme(plot.caption=element_text(hjust=0))
ggsave('forest_trends_over_time.png')

## Carbon Stock by Forest Age
carbonValsByAge <- rbindlist(lapply(carbonStockFiles,  function(x) readRDS(x)$csForest),idcol = 'windowid',fill=TRUE,use.names= TRUE)
carbonValsNonForest <- rbindlist(lapply(carbonStockFiles,  function(x) readRDS(x)$csNonForest),idcol = 'windowid',fill=TRUE)
                                        #carbonValsByAge[forest_ageV == Inf, forest_ageV := 40] ##Asign 40 to the forest that's oldest for the purposes of the graphs and regressions
nameCol <- which(names(carbonValsByAge) == 'pixelMapBiomasIndex')
if (length(nameCol)>1)
        set(carbonValsByAge,,nameCol[2:length(nameCol)],NULL) ##Take away duplicate pixelMapBiomasIndex

#Carbon as function of forest age graphs
ggplot(carbonValsByAge, aes(x=forest_ageV, y=carbonVal)) +
    stat_summary() +
    stat_smooth() + 
    theme_bw()+
    scale_x_continuous('Age')+
    scale_y_continuous('Tons of Carbon Per Hectare')+
    labs(caption = 'The carbon stock estimates are from Englund, O. et.al.\nAll pixels shown are conditional on being non-forest at some point since 1985.\nLoess trend line is in blue.')+
    theme(plot.caption=element_text(hjust=0))
ggsave('carbon_stock_forest_age.png')

## Regressions
carbonValsAll <- rbindlist(list(forest = carbonValsByAge, nonForest = carbonValsNonForest),
                           use.names = TRUE, fill = TRUE, idcol = 'landtyp')
carbonValsAll[,recodedLandUse.viterbiF := NULL]
carbonValsAll[is.na(forest_ageV), forest_ageV :=0]
carbonValsAll[,forest_age := forest_ageV] ##to simplify use of predict in the carbon simulation
carbonValsAll[, forestDum := landtyp == 'forest']
windowIDWithForest <- unique(carbonValsAll[landtyp == 'forest']$windowid)

##This is just meant to enable the running of the regressions (since it will drop rows with Inf
longForestAge <- 50
                                        #carbonRegVFE <- felm(carbonVal ~ I(forest_ageV * (forest_ageV>0))|factor(windowid)*factor(forestDum),carbonValsAll[windowid %in% windowIDWithForest])
carbonValsAll[is.infinite(forest_age), forest_age := longForestAge] ##This is just a dummy value so that it doesn't get eliminated from the regression
carbonRegV <- lm(carbonVal ~ I(forest_age * (forest_age>0 & forest_age!=longForestAge)) + forestDum + I(forest_age == longForestAge),carbonValsAll)
stargazer(carbonRegV,type = 'latex',
          covariate.labels = c('$\\tilde\\alpha$', '$\\tilde\\beta$','$\\tilde\\gamma$','$\\delta$'),
          dep.var.caption = "",
          dep.var.labels = 'Carbon Stock',
          notes = 'Data as described in text. Regression uses data from 2017.')

##Distribution of forest age by frequency and viterbi for 2017
landUse2017 <- rbindlist(lapply(carbonStockFiles,  function(x) readRDS(x)$landuse[year==2017,-1]),idcol = 'windowid',fill=TRUE)
datForPlt <- melt(landUse2017, measure = patterns("^forest_age"), id.vars = c('year','pixelMapBiomasIndex'),value.name = c('forest_age'))
datForPlt[is.infinite(forest_age),forest_age := longForestAge]
datForPlt[variable == 'forest_ageV',variable := 'Viterbi']
datForPlt[variable == 'forest_ageY',variable := 'Raw']

percentOfForestAtMax <- datForPlt[variable == 'Viterbi' & !is.na(forest_age),mean(forest_age == longForestAge)]

ggplot(datForPlt[forest_age>0 & variable %in% c('Viterbi','Raw')],aes(x = forest_age,color = variable)) +
    stat_ecdf()+
    theme_bw()+
    scale_x_continuous('Age')+
    scale_y_continuous('Share of Forest Under Age')+
    coord_cartesian(xlim = c(1,32))+
    scale_color_manual('Measurement',values = cbPalette)+
    ggtitle('Cumulative Distribution of Forest Age By Measurement Type',subtitle = '2017 Data')+
    labs(caption = 'Pixels that were continuously forest since 1985 are over 32 years old, but their age is otherwise unknown.')+
    theme(plot.caption=element_text(hjust=0))
ggsave('forest_age_dist.png')
rm(datForPlt)

##Carbon Stock Calculation
##Import landuse data for all years

carbonStockCalcFunc <- function(fileNm,carbonRegRes = carbonRegV){
    landUse <- readRDS(fileNm)$landuse
    landUse[,which(duplicated(names(landUse))) := NULL]

    
    setnames(landUse,names(landUse),str_replace(names(landUse),'_lag_0',''))

    landUseMelt <- melt(landUse, measure = patterns("^forestDum","^forest_age"), id.vars = c('year','pixelMapBiomasIndex'),value.name = c('forestDum','forest_age'))

    levels(landUseMelt$variable) <- c('Viterbi','Raw','Baseline','No Deforestation')
    landUseMelt[year <= 2020,mean(is.na(forestDum)),by = variable]
    landUseMelt[year > 2020,mean(forestDum,na.rm=TRUE),by = variable]
    
    landUseMelt[is.na(forest_age),forest_age := 0]
    landUseMelt[is.infinite(forest_age),forest_age := longForestAge]

    ##Carbon stock simulation
    landUseMelt[,carbonPred := predict(carbonRegRes,landUseMelt)]
    carbonStockByYear <- landUseMelt[,list(carbonPred = mean(carbonPred,na.rm=TRUE)), by = list(variable,year)]

    return(carbonStockByYear)
}

carbonStockResList <- pblapply(carbonStockFiles, carbonStockCalcFunc)
carbonStockByYear <- rbindlist(carbonStockResList)[,list(carbonPred = mean(carbonPred) * hectaresAtlanticForest),
                                                   by = list(variable,year)]

##Compute 2020 carbon stock and convert to potential CO2 emissions

carbonStock2020 <- carbonStockByYear[year == 2020 & variable %in% c('Viterbi','Raw'),list(Measurement = variable, `Carbon Stock (Bn Tons)` = round(carbonPred/1e9,2), `Social Value (in Bn Dollars)` = round(carbonPred*carbonToDollars2020/1e9,2))]
knitr::kable(carbonStock2020,'latex')

carbonStock2040 <- carbonStockByYear[year == 2040 & variable %in% c('Baseline','No Deforestation'),list(Measurement = variable, `Carbon Stock (Bn Tons)` = round(carbonPred/1e9,2), `Social Value (in Bn Dollars)` = round(carbonPred*carbonToDollars2040/1e9,2))]
knitr::kable(carbonStock2040,'latex')
