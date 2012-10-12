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

fullMat <- cbind(rawDatMatList$zhu[intersectFeatures, ],
                 rawDatMatList$hou[intersectFeatures, ],
                 rawDatMatList$dir[intersectFeatures, ],
                 rawDatMatList$lusc[intersectFeatures, ])

studyIndicator <- c(rep('zhu', ncol(rawDatMatList$zhu)),
                    rep('hou', ncol(rawDatMatList$hou)),
                    rep('dir', ncol(rawDatMatList$dir)),
                    rep('lusc', ncol(rawDatMatList$lusc)))

svdObj <- svd(fullMat)

rawDF <- data.frame(svdObj$v[ , 1:2])
colnames(rawDF) <- c('PrinComp1', 'PrinComp2')

rawPcPlot <- ggplot(rawDF, aes(PrinComp1, PrinComp2)) +
  geom_point(aes(colour = factor(studyIndicator),
                 shape = factor(studyIndicator),
                 size = 20)) +
                   scale_size(guide = 'none')
                    





