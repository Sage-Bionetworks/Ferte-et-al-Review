## figureFour.R

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washington
## erich.huang@sagebase.org

## Heatmap visualiations of highly-weighted features in the models

## REQUIRED LIBRARIES
require(synapseClient)
require(ggplot2)
ggHeatEnt <- loadEntity('syn274063')
attach(ggHeatEnt)

## PULL DOWN THE predictiveModels.R WORKSPACE
predModEnt <- loadEntity('syn1447949')

# Unpack the workspace from the entity
objNames <- names(predModEnt$objects)
for(i in 1:length(predModEnt$objects)){
  varName <- objNames[i]
  targetObj <- predModEnt$objects[[i]]
  assign(varName, targetObj)
}