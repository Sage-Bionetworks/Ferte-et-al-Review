## figureOnePanelSix.R

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

## FIFTH, PLOTTING FRMA NORMALIZED DATA
zhuFrmaEnt <- loadEntity('syn1437066')
houFrmaEnt <- loadEntity('syn1437182')
dirFrmaEnt <- loadEntity('syn1437194')
luscFrmaEnt <- loadEntity('syn1437182')

frmaEntList <- list('zhu' = zhuFrmaEnt,
                     'hou' = houFrmaEnt,
                     'dir' = dirFrmaEnt,
                     'lusc' = luscFrmaEnt)

frmaDatList <- lapply(frmaEntList, function(x){
  exprs <- exprs(x$objects[[1]])
})

studyIndicator <- c(rep('zhu', ncol(frmaDatList$zhu)), ##
                    rep('hou', ncol(frmaDatList$hou)), ##
                    rep('dir', ncol(frmaDatList$dir)), ##
                    rep('lusc', ncol(frmaDatList$lusc))) ##

commonFeatures <- intersectFeatures(frmaDatList) ##

intFrmaDatList <- lapply(frmaDatList, function(x){x[commonFeatures, ]})
fullFrmaMat <- Reduce(cbind, intFrmaDatList)

frmaPcPlot <- generatePcPlot(fullFrmaMat) +
  opts(title = 'Separate fRMA Normalization by Prin. Comp.\n')
frmaPcPlot