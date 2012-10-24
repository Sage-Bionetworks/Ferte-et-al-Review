## figureOnePanelFour.R

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

## FOURTH, PLOTTING GCRMA NORMALIZED DATA
zhuGcrmaEnt <- loadEntity('syn1437007')
houGcrmaEnt <- loadEntity('syn1437176')
dirGcrmaEnt <- loadEntity('syn1437188')
luscGcrmaEnt <- loadEntity('syn1445137')

gcrmaEntList <- list('zhu' = zhuGcrmaEnt,
                     'hou' = houGcrmaEnt,
                     'dir' = dirGcrmaEnt,
                     'lusc' = luscGcrmaEnt)

gcrmaDatList <- lapply(gcrmaEntList, function(x){
  exprs <- exprs(x$objects[[1]])
})

studyIndicator <- c(rep('zhu', ncol(gcrmaDatList$zhu)), ##
                    rep('hou', ncol(gcrmaDatList$hou)), ##
                    rep('dir', ncol(gcrmaDatList$dir)), ##
                    rep('lusc', ncol(gcrmaDatList$lusc))) ##

commonFeatures <- intersectFeatures(gcrmaDatList) ##

intGcrmaDatList <- lapply(gcrmaDatList, function(x){x[commonFeatures, ]})
fullGcrmaMat <- Reduce(cbind, intGcrmaDatList)

gcrmaPcPlot <- generatePcPlot(fullGcrmaMat) +
  opts(title = 'Separate GCRMA Normalization by Prin. Comp.\n')
gcrmaPcPlot