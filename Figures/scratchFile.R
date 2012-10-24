## Scratch temp file

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washington
## erich.huang@sagebase.org

## FIRST: GENERATING UN-NORMALIZED, UN-BACKGROUND CORRECTED DATA
## We'll use the 'rma()' function with the normalize and background arguments
## set to false.

## PULL IN THE RAW DATA FROM SYNAPSE
zhuRawEnt <- loadEntity('syn1439020')
houRawEnt <- loadEntity('syn1422295')
dirRawEnt <- loadEntity('syn1422422')
luscRawEnt <- loadEntity('syn1426948')

rawEntities <- list(zhuRawEnt, houRawEnt, dirRawEnt, luscRawEnt)
names(rawEntities) <- c('zhu', 'hou', 'dir', 'lusc')

celNamesList <- lapply(rawEntities, getCelNames)

rawDatList <- lapply(celNamesList, getRawDat)

rawDatMatList <- lapply(rawDatList, exprs)

# Since these platforms are on two different Affymetrix platforms, we use
# the intersection of the Affymetrix probesets for comparative purposes

commonFeatures <- intersectFeatures(rawDatMatList)

fullRawMat <- cbind(rawDatMatList$zhu[intersectFeatures, ],
                 rawDatMatList$hou[intersectFeatures, ],
                 rawDatMatList$dir[intersectFeatures, ],
                 rawDatMatList$lusc[intersectFeatures, ])

studyIndicator <- c(rep('zhu', ncol(rawDatMatList$zhu)),
                    rep('hou', ncol(rawDatMatList$hou)),
                    rep('dir', ncol(rawDatMatList$dir)),
                    rep('lusc', ncol(rawDatMatList$lusc)))

# Plot the raw expression data
rawPcPlot <- generatePcPlot(fullRawMat) + 
  opts(title = 'Raw without Normalization or Background Correction by Prin. Comp.\n')
rawPcPlot

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

intRmaDatList <- lapply(rmaDatList, function(x){x[intersectFeatures, ]})
fullRmaMat <- Reduce(cbind, intRmaDatList)

rmaPcPlot <- generatePcPlot(fullRmaMat) + 
  opts(title = 'Separate RMA Normalization by Prin. Comp.\n')
rmaPcPlot

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

intGcrmaDatList <- lapply(gcrmaDatList, function(x){x[intersectFeatures, ]})
fullGcrmaMat <- Reduce(cbind, intGcrmaDatList)

gcrmaPcPlot <- generatePcPlot(fullGcrmaMat) +
  opts(title = 'Separate GCRMA Normalization by Prin. Comp.\n')
gcrmaPcPlot

## FIFTH, PLOTTING DCHIP NORMALIZED DATA
zhuDchipEnt <- loadEntity('syn1437063')
houDchipEnt <- loadEntity('syn1437180')
dirDchipEnt <- loadEntity('syn1437192')
luscDchipEnt <- loadEntity('syn1437118')

dchipEntList <- list('zhu' = zhuDchipEnt,
                     'hou' = houDchipEnt,
                     'dir' = dirDchipEnt,
                     'lusc' = luscDchipEnt)

dchipDatList <- lapply(dchipEntLIst, function(x){
  exprs <- exprs(x$objects[[1]])
})

intDchipDatList <- lapply(dchipDatList, function(x){x[intersectFeatures, ]})
fullDchipMat <- Reduce(cbind, intDchipDatList)

dchipPcPlot <- generatePcPlot(fullDchipMat) +
  opts(title = 'Separate dCHIP Normalization by Prin. Comp.\n')
dchipPcPlot

## PUT ALL THE PLOTS TOGETHER
# Source in a multiplot function
multiplotEnt <- loadEntity('syn274067')
attach(multiplotEnt)
fullFigureOne <- multiplot(rawPcPlot, 
                           mas5PcPlot, 
                           rmaPcPlot, 
                           gcrmaPcPlot,
                           cols = 2)
