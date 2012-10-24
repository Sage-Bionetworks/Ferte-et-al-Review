## figureOnePanelFive.R

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

## FIFTH, PLOTTING DCHIP NORMALIZED DATA
zhuDchipEnt <- loadEntity('syn1437063')
houDchipEnt <- loadEntity('syn1437180')
dirDchipEnt <- loadEntity('syn1437192')
luscDchipEnt <- loadEntity('syn1437118')

dchipEntList <- list('zhu' = zhuDchipEnt,
                     'hou' = houDchipEnt,
                     'dir' = dirDchipEnt,
                     'lusc' = luscDchipEnt)

dchipDatList <- lapply(dchipEntList, function(x){
  exprs <- exprs(x$objects[[1]])
})

studyIndicator <- c(rep('zhu', ncol(dchipDatList$zhu)), ##
                    rep('hou', ncol(dchipDatList$hou)), ##
                    rep('dir', ncol(dchipDatList$dir)), ##
                    rep('lusc', ncol(dchipDatList$lusc))) ##

commonFeatures <- intersectFeatures(dchipDatList) ##

intDchipDatList <- lapply(dchipDatList, function(x){x[commonFeatures, ]})
fullDchipMat <- Reduce(cbind, intDchipDatList)

dchipPcPlot <- generatePcPlot(fullDchipMat) +
  opts(title = 'Separate dCHIP Normalization by Prin. Comp.\n')
dchipPcPlot