## figureThree.R

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washington
## erich.huang@sagebase.org

## Visualizing the results of clinico-genomic models of 3-year overall survival in lung cancer

## REQUIRED LIBRARIES
require(synapseClient)
require(ggplot2)
require(survival)

## PULL DOWN THE predictiveModels.R WORKSPACE
predModEnt <- loadEntity('syn1447949')

# Unpack the workspace from the entity
objNames <- names(predModEnt$objects)
for(i in 1:length(predModEnt$objects)){
  varName <- objNames[i]
  targetObj <- predModEnt$objects[[i]]
  assign(varName, targetObj)
}

foo <- survfit(Surv(zhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhuClin$VITAL_STATUS) ~ riskClin)