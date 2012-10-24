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
entityIdList <- as.list(figureOneQuery[ , 2])

figOneEntList <- lapply(entityIdList, loadEntity)
names(figOneEntList) <- figureOneQuery[ , 1]

figOneObjList <- lapply(figOneEntList, function(x){x$objects[[1]]})
studyIndicator <- figOneObjList$figureOneStudyIndicator

figOneObjList$figureOnePanelOneBinary
figOneObjList$figureOnePanelTwoBinary
figOneObjList$figureOnePanelThreeBinary
figOneObjList$figureOnePanelFourBinary
figOneObjList$figureOnePanelFiveBinary
figOneObjList$figureOnePanelSixBinary

multiplot(figOneObjList, cols = 3)



