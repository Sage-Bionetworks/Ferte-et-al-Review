## figureOnePanelThree.R

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washington
## erich.huang@sagebase.org

## REQUIRED LIBRARIES
require(Biobase)
require(affy)
require(snm)
require(ggplot2)
require(synapseClient)

figFuncEnt <- loadEntity('syn1446521')
attach(figFuncEnt)

## THIRD, PLOTTING RMA NORMALIZED DATA
## PULL IN RMA DATA FROM SYNAPSE
zhuRmaEnt <- loadEntity('syn1436971')
houRmaEnt <- loadEntity('syn1437174')
dirRmaEnt <- loadEntity('syn1440819')
luscRmaEnt <- loadEntity('syn1437109')

rmaEntList <- list('zhu' = zhuRmaEnt,
                   'hou' = houRmaEnt,
                   'dir' = dirRmaEnt,
                   'lusc' = luscRmaEnt)

rmaDatList <- lapply(rmaEntList, function(x){
  exprs <- exprs(x$objects[[1]])
})

studyIndicator <- c(rep('zhu', ncol(rmaDatList$zhu)), ##
                    rep('hou', ncol(rmaDatList$hou)), ##
                    rep('dir', ncol(rmaDatList$dir)), ##
                    rep('lusc', ncol(rmaDatList$lusc))) ##

commonFeatures <- intersectFeatures(rmaDatList) ##

intRmaDatList <- lapply(rmaDatList, function(x){x[commonFeatures, ]}) ##
fullRmaMat <- Reduce(cbind, intRmaDatList)

rmaPcPlot <- generatePcPlot(fullRmaMat) + 
  opts(title = 'Separate RMA Normalization by Prin. Comp.\n')
rmaPcPlot