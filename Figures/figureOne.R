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

## FIRST: EXTRACTING RAW DATA FROM THE CEL FILES
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

fullMat <- cbind(rawDatList$zhu, 
                 rawDatList$hou, 
                 rawDatList$dir, 
                 rawDatList$lusc)



