## figureOnePanelSeven.R

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washington
## erich.huang@sagebase.org

## REQUIRED LIBRARIES
require(synapseClient)

## FIND THE SYNAPSE ENTITIES ON HGU133A PLATFORM
projectQuery <- synapseQuery('SELECT id, name FROM entity WHERE entity.parentId == "syn87682"')
hgu133aEntInd <- grep('HGU133a', projectQuery[ , 1])
hgu133aEntIds <- as.list(projectQuery[hgu133aEnt])

hgu133aEntList <- lapply()


