## displayFigureOnePanels.R

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washington
## erich.huang@sagebase.org

## REQUIRE LIBRARIES
require(synapseClient)
multiplotEnt <- loadEntity('syn274067')
attach(multiplotEnt)

## QUERY FIGURE ONE FOLDER
figureOneQuery <- synapseQuery('SELECT id, name FROM entity WHERE entity.parentId == "syn1446518"')

## LAPPLY THROUGH THE QUERY DATAFRAME AND LOAD ENTITIES
entityIdList <- as.list(figureOneQuery[2:8, 2])

figOneEntList <- lapply(entityIdList, loadEntity)
names(figOneEntList) <- figureOneQuery[2:8, 1]

figOneObjList <- lapply(figOneEntList, function(x){x$objects[[1]]})
studyIndicator <- figOneObjList$figureOneStudyIndicator

# figOneObjList$figureOnePanelOneBinary
# figOneObjList$figureOnePanelTwoBinary
# figOneObjList$figureOnePanelThreeBinary
# figOneObjList$figureOnePanelFourBinary
# figOneObjList$figureOnePanelFiveBinary
# figOneObjList$figureOnePanelSixBinary

## DISPLAY ALL THE PLOTS IN A 3X2 PANEL
multiplot(plotlist = figOneObjList[1:6], cols = 3)

## FOR THE SECOND PART OF FIGURE ONE, ALL THE BINARIES ARE STORED IN ONE LIST OBJECT
figOneHgu133aEnt <- loadEntity('syn1446680')
figOneHgu133aObjs <- figOneHgu133aEnt$objects$figureOneHgu133aPanelList

studyIndicator <- figOneHgu133aObjs$studyInd

## DISPLAY THESE PLOTS IN A 3X2 PANEL
multiplot(plotlist = figOneHgu133aObjs[1:6], cols = 3)






