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

getCelNames <- function(x){
  fileNames <- list.celfiles(path = x$cacheDir, full.names = TRUE)
}

rawEntities <- list(zhuRawEnt, houRawEnt, dirRawEnt, luscRawEnt)
names(rawEntities) <- c('zhu', 'hou', 'dir', 'lusc')

celNamesList <- lapply(rawEntities, getCelNames)

getRawDat <- function(x){
  affyBatchObj <- ReadAffy(filenames = x)
  rawDat <- rma(affyBatchObj, normalize = FALSE, background = FALSE)
}

rawDatList <- lapply(celNamesList, getRawDat)

rawDatMatList <- lapply(rawDatList, exprs)

featureList <- lapply(rawDatMatList, rownames)

intersectFeatures <- Reduce(intersect, featureList)

fullRawMat <- cbind(rawDatMatList$zhu[intersectFeatures, ],
                 rawDatMatList$hou[intersectFeatures, ],
                 rawDatMatList$dir[intersectFeatures, ],
                 rawDatMatList$lusc[intersectFeatures, ])

studyIndicator <- c(rep('zhu', ncol(rawDatMatList$zhu)),
                    rep('hou', ncol(rawDatMatList$hou)),
                    rep('dir', ncol(rawDatMatList$dir)),
                    rep('lusc', ncol(rawDatMatList$lusc)))

rawSvdObj <- fast.svd(fullRawMat)

rawDF <- data.frame(rawSvdObj$v[ , 1:2])
colnames(rawDF) <- c('PrinComp1', 'PrinComp2')

rawPcPlot <- ggplot(rawDF, aes(PrinComp1, PrinComp2)) +
  geom_point(aes(colour = factor(studyIndicator),
                 size = 20)) +
                   scale_size(guide = 'none')

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

mas5SvdObj <- fast.svd(fullMas5Mat)
mas5DF <- data.frame(mas5SvdObj$v[ , 1:2])
colnames(mas5DF) <- c('PrinComp1', 'PrinComp2')

mas5PcPlot <- ggplot(mas5DF, aes(PrinComp1, PrinComp2)) +
  geom_point(aes(colour = factor(studyIndicator),
                 size = 20)) +
                   scale_size(guide = 'none')

## THIRD, PLOTTING RMA NORMALIZED DATA
## PULL IN RMA DATA FROM SYNAPSE
zhuRmaEnt <- loadEntity('syn1436971')
houRmaEnt <- loadEntity('syn1437174')
dirRmaEnt <- loadEntity('syn1437186')
luscRmaEnt <- loadEntity('syn1437109')

rmaEntList <- list('zhu' = zhuRmaEnt,
                   'hou' = houRmaEnt,
                   'dir' = dirRmaEnt,
                   'lusc' = luscRmaEnt)

rmaDatList <- lapply(rmaEntList, function(x){
  exprs <- exprs(x$objects[[1]])
})