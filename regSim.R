library(Rsolnp)
library(data.table)
library(grid)
library(parallel)
library(plyr)

##Run MC simulations using regression
source("hmm_functions.R")
source("hmm_parameters.R")

cl <- makeCluster(8)

nCells <- 100
nPixelPerCell <- 1000
alpha <- -3
beta <- 1
nSim <- 100

set.seed(321321)

xVec <- rnorm(nCells)
trueDeforestPr <- exp(alpha + beta * xVec)/(1+exp(alpha + beta * xVec))
tMax <- length(params0$P_list)+1


##Deforestation only varies in the 3rd period
paramList <- list()
for (i in 1:length(xVec)){
    paramList[[i]] <- params0
    paramList[[i]]$P_list[[3]][1,2] <- trueDeforestPr[i]
    paramList[[i]]$P_list[[3]][1,1] <- 1 - trueDeforestPr[i]
}

##Up until here is just simulating the dataset, so it is the same across MC iterations
regResList <- list()
for (s in 1:nSim){
    message('Starting iteration ',s)
    datDraw <- mclapply(paramList,
                        function(x) replicate(nPixelPerCell, simulate_hmm(x), simplify=FALSE),mc.cores = 8)

    params0_hat <- mclapply(datDraw,
                            function(x) get_expectation_maximization_estimates(x, params0, max_iter=20, epsilon=0.001),mc.cores=8)

    viterbiList <- lapply(1:length(xVec),
                          function(i) apply_viterbi_path_in_parallel(datDraw[[i]], params_hat=params0_hat[[i]], max_cores=8))

    ##Reshape data for estimation
    dataSetList <- list()
    for (i in 1:length(xVec)){
        dataSetList[[i]] <- list()
        for (j in 1:nPixelPerCell){
            dataSetList[[i]][[j]] <- data.table(t = 1:tMax, viterbiVal = viterbiList[[i]][[j]],
                                                yVal = datDraw[[i]][[j]]$y,
                                                trueSVal = datDraw[[i]][[j]]$x)
        }
        dataSetList[[i]] <- rbindlist(dataSetList[[i]],idcol='pixelIDTemp')
    }

    dataSet <- rbindlist(dataSetList,idcol = 'cellID')
    dataSet[,pixelID := (cellID-1)*nPixelPerCell + pixelIDTemp]
    dataSet[,pixelIDTemp := NULL]

    dataSetMelt <- melt(dataSet,id.vars = c('cellID','pixelID','t'), measure.vars = c('viterbiVal','yVal','trueSVal'))

    dataSetMelt[, valueLag := shift(value,1),by = list(variable,cellID,pixelID)]

    regressionData <- dataSetMelt[t == 4 & valueLag == 1, list(deforestRate = mean(value==2)),
                                  by = list(cellID,variable)]

    regressionData <- merge(regressionData,
                            data.table(cellID = 1:length(xVec),xVal = xVec))

    ##Run regression
    regResList[[s]] <- dlply(regressionData,
                    'variable',
                    function(dat) glm(deforestRate ~ xVal, family = 'binomial',data =dat))

}
