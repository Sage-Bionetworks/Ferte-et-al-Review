## Andrew Trister
## Sage Bionetworks


## Raw Data PC plots - modified from Erich Huang's code
## FIRST: GENERATING UN-NORMALIZED, UN-BACKGROUND CORRECTED DATA
## We'll use the 'rma()' function with the normalize and background arguments
## set to false.

## REQUIRED LIBRARIES
require(Biobase)
require(affy)
require(snm)
require(ggplot2)
require(synapseClient)


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

fullRawMat <- cbind(rawDatMatList$zhu[commonFeatures, ],
                    rawDatMatList$hou[commonFeatures, ],
                    rawDatMatList$dir[commonFeatures, ],
                    rawDatMatList$lusc[commonFeatures, ])

studyIndicator <- c(rep('zhu', ncol(rawDatMatList$zhu)),
                    rep('hou', ncol(rawDatMatList$hou)),
                    rep('dir', ncol(rawDatMatList$dir)),
                    rep('lusc', ncol(rawDatMatList$lusc)))

#Plot the principal components

svdObj <- fast.svd(fullRawMat)
plot(svdObj$v[,1],svdObj$v[,2],
     col=c("royal blue","orange","aquamarine4","brown2")[as.factor(studyIndicator)],
     bg=c("royal blue","orange","aquamarine4","brown2")[as.factor(studyIndicator)],
     pch=c(8,22,23,24)[as.factor(studyIndicator)],
     cex=0.8)


