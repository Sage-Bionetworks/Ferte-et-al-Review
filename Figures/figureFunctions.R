## figureFunctions.R

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washington
## erich.huang@sagebase.org

## R functions for generating figures.

## REQUIRED LIBRARIES
require(Biobase)
require(affy)
require(snm)
require(ggplot2)
require(synapseClient)

## synapseLogin()


## DEFINE FUNCTIONS
## DEFINE getCelNames() FUNCTION
getCelNames <- function(x){
  fileNames <- list.celfiles(path = x$cacheDir, full.names = TRUE)
}

## DEFINE getRawDat() FUNCTION
# Merely using the rma() function to extract summarized but not normalized or
# background-adjusted expression data
getRawDat <- function(x){
  affyBatchObj <- ReadAffy(filenames = x)
  rawDat <- rma(affyBatchObj, normalize = FALSE, background = FALSE)
}

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

## DEFINE intersectFeatures() FUNCTION
## Take the entity list, get the feature names for the disparate Affy platforms and return
## the intersection of the feature names across platforms
intersectFeatures <- function(datMatList){
  featureList <- lapply(datMatList, rownames)
  intersectNames <- Reduce(intersect, featureList)
  return(intersectNames)
}
