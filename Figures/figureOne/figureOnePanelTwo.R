## figureOnePanelTwo.R

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

## SECOND, PLOTTING MAS5 NORMALIZED DATA
## PULL IN MAS5 DATA FROM SYNAPSE
zhuMas5Ent <- loadEntity('syn1437053')
houMas5Ent <- loadEntity('syn1437178')
dirMas5Ent <- loadEntity('syn1437190')
luscMas5Ent <- loadEntity('syn1437114')

mas5EntList <- list('zhu' = zhuMas5Ent,
                    'hou' = houMas5Ent,
                    'dir' = dirMas5Ent,
                    'lusc' = luscMas5Ent)

mas5DatList <- lapply(mas5EntList, function(x){
  exprs <- exprs(x$objects[[1]])
})

intMas5DatList <- lapply(mas5DatList, function(x){x[intersectFeatures, ]})
fullMas5Mat <- Reduce(cbind, intMas5DatList)

mas5PcPlot <- generatePcPlot(fullMas5Mat) +
  opts(title = 'Separate MAS5 Normalization by Prin. Comp.\n')
mas5PcPlot