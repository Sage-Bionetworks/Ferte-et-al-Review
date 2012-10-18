## figureOne.R

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washington
## erich.huang@sagebase.org

## R script for generating figure one of the paper.

## REQUIRED LIBRARIES
require(Biobase)
require(affy)
require(snm)
require(ggplot2)
require(synapseClient)

## synapseLogin()

## FIRST: GENERATING UN-NORMALIZED, UN-BACKGROUND CORRECTED DATA
## We'll use the 'rma()' function with the normalize and background arguments
## set to false.

## PULL IN THE RAW DATA FROM SYNAPSE
zhuRawEnt <- loadEntity('syn1439020')
houRawEnt <- loadEntity('syn1422295')
dirRawEnt <- loadEntity('syn1422422')
luscRawEnt <- loadEntity('syn1426948')

## DEFINE getCelNames() FUNCTION
getCelNames <- function(x){
  fileNames <- list.celfiles(path = x$cacheDir, full.names = TRUE)
}

rawEntities <- list(zhuRawEnt, houRawEnt, dirRawEnt, luscRawEnt)
names(rawEntities) <- c('zhu', 'hou', 'dir', 'lusc')

celNamesList <- lapply(rawEntities, getCelNames)

## DEFINE getRawDat() FUNCTION
# Merely using the rma() function to extract summarized but not normalized or
# background-adjusted expression data
getRawDat <- function(x){
  affyBatchObj <- ReadAffy(filenames = x)
  rawDat <- rma(affyBatchObj, normalize = FALSE, background = FALSE)
} 

rawDatList <- lapply(celNamesList, getRawDat)

rawDatMatList <- lapply(rawDatList, exprs)

featureList <- lapply(rawDatMatList, rownames)

# Since these platforms are on two different Affymetrix platforms, we use
# the intersection of the Affymetrix probesets for comparative purposes

intersectFeatures <- Reduce(intersect, featureList)

fullRawMat <- cbind(rawDatMatList$zhu[intersectFeatures, ],
                 rawDatMatList$hou[intersectFeatures, ],
                 rawDatMatList$dir[intersectFeatures, ],
                 rawDatMatList$lusc[intersectFeatures, ])

studyIndicator <- c(rep('zhu', ncol(rawDatMatList$zhu)),
                    rep('hou', ncol(rawDatMatList$hou)),
                    rep('dir', ncol(rawDatMatList$dir)),
                    rep('lusc', ncol(rawDatMatList$lusc)))

## DEFINE generatePcPlot() FUNCTION
# Create an SVD object, make a dataframe out of the first two factors
# and plot it using ggplot
generatePcPlot <- function(fullMatrix){
  require(ggplot2)
  ## Make a ggplot-friendly dataframe from the singular value decomposition
  svdObj <- fast.svd(fullMatrix)
  pcDF <- data.frame(svdObj$v[ , 1:2])
  colnames(pcDF) <- c('PrinComp1', 'PrinComp2')
  ## Plot it
  pcPlot <- ggplot(pcDF, aes(PrinComp1, PrinComp2)) +
    geom_point(aes(colour = factor(studyIndicator), size = 20)) +
    scale_size(guide = 'none')
}

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
  opts(title = 'Seperate RMA Normalization by Prin. Comp.\n')
rmaPcPlot

## FOURTH, PLOTTING GCRMA NORMALIZED DATA
zhuGcrmaEnt <- loadEntity('syn1437007')
houGcrmaEnt <- loadEntity('syn1437176')
dirGcrmaEnt <- loadEntity('syn1437188')
luscGcrmaEnt <- loadEntity('syn1437111')

gcrmaEntList <- list('zhu' = zhuGcrmaEnt,
                     'hou' = houGcrmaEnt,
                     'dir' = dirGcrmaEnt,
                     'lusc' = luscGcrmaEnt)

gcrmaDatList <- lapply(gcrmaEntList, function(x){
  exprs <- exprs(x$objects[[1]])
})

## PUT ALL THE PLOTS TOGETHER
# Source in a multiplot function
multiplotEnt <- loadEntity('syn274067')
attach(multiplotEnt)
fullFigureOne <- multiplot(rawPcPlot, mas5PcPlot, rmaPcPlot, cols = 2)
