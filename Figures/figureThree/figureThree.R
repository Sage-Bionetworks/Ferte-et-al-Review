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
ggkmEnt <- loadEntity('syn299147')
attach(ggkmEnt)
require(plyr)

## PULL DOWN THE predictiveModels.R WORKSPACE
predModEnt <- loadEntity('syn1447949')

# Unpack the workspace from the entity
objNames <- names(predModEnt$objects)
for(i in 1:length(predModEnt$objects)){
  varName <- objNames[i]
  targetObj <- predModEnt$objects[[i]]
  assign(varName, targetObj)
}

## CREATE SURVFIT OBJECTS
# Clinical only
clinSurvFit <- survfit(Surv(zhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH, 
                            zhuClin$VITAL_STATUS) ~ riskClin)

clinKmPlot <- ggkm(clinSurvFit,
                   returns = TRUE,
                   table = FALSE,
                   ylabs = 'Survival Probability\n',
                   xlabs = 'Time (Months)',
                   ystratalabs = c('High-Risk', 'Low-Risk'),
                   timeby = 50,
                   main = 'Clinical Logit Model\n')

clinKmPlot <- clinKmPlot + 
  geom_vline(xintercept = 36, linetype = 'longdash', colour = 'red', size = 0.125) +
  opts(legend.position = c(0.9, 0.9))

# Elasticnet
enetSurvFit <- survfit(Surv(zhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH,
                            zhuClin$VITAL_STATUS) ~ riskEnet)

enetKmPlot <- ggkm(enetSurvFit,
                   returns = TRUE,
                   table = FALSE,
                   ylabs = 'Survival Probability\n',
                   ystratalabs = c('High-Risk', 'Low-Risk'),
                   timeby = 50,
                   main = 'Elasticnet Model\n')

enetKmPlot <- enetKmPlot +
  geom_vline(xintercept = 36, linetype = 'longdash', colour = 'red', size = 0.125)+
  opts(legend.position = c(0.9, 0.9))

# Principle Component Regression
pcrSurvFit <- survfit(Surv(zhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH,
                            zhuClin$VITAL_STATUS) ~ riskPcr)

pcrKmPlot <- ggkm(pcrSurvFit,
                   returns = TRUE,
                   table = FALSE,
                   ylabs = 'Survival Probability\n',
                   ystratalabs = c('High-Risk', 'Low-Risk'),
                   timeby = 50,
                   main = 'Principle Component Regression Model\n')

pcrKmPlot <- pcrKmPlot +
  geom_vline(xintercept = 36, linetype = 'longdash', colour = 'red', size = 0.125)+
  opts(legend.position = c(0.9, 0.9))

# Partial Least Squares Regression
plsSurvFit <- survfit(Surv(zhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH,
                           zhuClin$VITAL_STATUS) ~ riskPls)

plsKmPlot <- ggkm(plsSurvFit,
                  returns = TRUE,
                  table = FALSE,
                  ylabs = 'Survival Probability\n',
                  ystratalabs = c('High-Risk', 'Low-Risk'),
                  timeby = 50,
                  main = 'Principle Component Regression Model\n')

plsKmPlot <- plsKmPlot +
  geom_vline(xintercept = 36, linetype = 'longdash', colour = 'red', size = 0.125)+
  opts(legend.position = c(0.9, 0.9))

# Random Forest
rfSurvFit <- survfit(Surv(zhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH,
                           zhuClin$VITAL_STATUS) ~ riskRF)

rfKmPlot <- ggkm(rfSurvFit,
                  returns = TRUE,
                  table = FALSE,
                  ylabs = 'Survival Probability\n',
                  ystratalabs = c('High-Risk', 'Low-Risk'),
                  timeby = 50,
                  main = 'Principle Component Regression Model\n')

rfKmPlot <- rfKmPlot +
  geom_vline(xintercept = 36, linetype = 'longdash', colour = 'red', size = 0.125)+
  opts(legend.position = c(0.9, 0.9))