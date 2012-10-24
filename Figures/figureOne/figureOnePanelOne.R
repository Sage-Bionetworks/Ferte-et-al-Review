## figureOnePanelOne.R

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washington
## erich.huang@sagebase.org

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