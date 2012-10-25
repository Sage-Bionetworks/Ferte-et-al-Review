## figureOnePanelSeven.R

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washington
## erich.huang@sagebase.org

## REQUIRED LIBRARIES
require(synapseClient)
figFuncEnt <- loadEntity('syn1446521')
attach(figFuncEnt)

## FIND THE SYNAPSE ENTITIES ON HGU133A PLATFORM
projectQuery <- synapseQuery('SELECT id, name FROM entity WHERE entity.parentId == "syn87682"')
projectQuery[ , 1] <- gsub(' ', '_', projectQuery[ , 1]) ## remove spaces in entity names

hgu133aEntInd <- grep('HGU133a', projectQuery[ , 1])
hgu133aEntIds <- as.list(projectQuery[hgu133aEntInd, 2])
hgu133aEntList <- lapply(hgu133aEntIds, loadEntity)

names(hgu133aEntList) <- projectQuery[hgu133aEntInd, 1]

## GET CLINICAL METADATA FOR THESE DATA
dirClinEnt <- loadEntity('syn1438222')
zhuClinEnt <- loadEntity('syn1438225')

## obtain the sample identifiers
dirSamps <- dirClinEnt$objects$DirClinF[ , 1]
zhuSamps <- zhuClinEnt$objects$ZhuClinF[ , 1]

zhuSamps <- paste(zhuSamps, '.CEL', sep = '')
studyIndicator <- c(rep('dir', length(dirSamps)), rep('zhu', length(zhuSamps)))

## RMA
rmaDat <- exprs(hgu133aEntList$HGU133a_rma$objects$HGU133A_rma)
rmaDat <- cbind(rmaDat[ , dirSamps], rmaDat[ , zhuSamps])

hgu133aRmaPcPlot <- generatePcPlot(rmaDat) +
  opts(title = 'RMA Normalized Together by Prin. Comp.\n')

## GCRMA
gcrmaDat <- exprs(hgu133aEntList$HGU133a_gcrma$objects$HGU133A_gcrma)
gcrmaDat <- cbind(gcrmaDat[ , dirSamps], gcrmaDat[ , zhuSamps])

hgu133aGcrmaPcPlot <- generatePcPlot(gcrmaDat) +
  opts(title = 'GCRMA Normalized Together by Prin. Comp.\n')

## MAS5
mas5Dat <- exprs(hgu133aEntList$HGU133a_MAS5$objects$HGU133A_MAS5)
mas5Dat <- cbind(mas5Dat[ , dirSamps], mas5Dat[ , zhuSamps])

hgu133aMas5PcPlot <- generatePcPlot(mas5Dat) +
  opts(title = 'MAS5 Normalized Together by Prin. Comp.\n')

## dCHIP
dchipDat <- exprs(hgu133aEntList$HGU133a_dCHIP$objects$HGU133A_dCHIP)
dchipDat <- cbind(dchipDat[ , dirSamps], dchipDat[ , zhuSamps])

hgu133aDchipPcPlot <- generatePcPlot(dchipDat) +
  opts(title = 'MAS5 Normalized Together by Prin. Comp.\n')

## BARCODE
barDat <- hgu133aEntList$HGU133a_barcode$objects$HGU133A_barcode
barDat <- cbind(barDat[ , dirSamps], barDat[ , zhuSamps])

hgu133aBarPcPlot <- generatePcPlot(barDat) +
  opts(title = 'Barcode Normalized by Prin. Comp.\n')

## FRMA
frmaDat <- exprs(hgu133aEntList$HGU133a_frma$objects$HGU133A_frma)
frmaDat <- cbind(frmaDat[ , dirSamps], frmaDat[ , zhuSamps])

hgu133aFrmaPcPlot <- generatePcPlot(frmaDat) +
  opts(title = 'fRMA Normalized by Prin. Comp.\n')
