## figureTwo.R

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washington
## erich.huang@sagebase.org

## Visualizing the results of clinico-genomic models of 3-year overall survival in lung cancer

## REQUIRED LIBRARIES
require(synapseClient)
require(ggplot2)

## ROC CURVES USING THE pROC OBJECTS
# Define a function 'makeRocDF' to make a dataframe suitable for ggplot2
makeRocDF <- function(rocObject){
  rocDF <- data.frame(rep(deparse(substitute(rocObject)), length(rocObject$specificities)),
                      1 - rocObject$specificities, rocObject$sensitivities)
  names(rocDF) <- c('Study', 'X', 'Y')
  return(rocDF)
}

foo <- makeRocDF(rocClin)

clinDF <- data.frame(rep('Clinical', length(rocClin$specificities)),
                         1 - rocClin$specificities, rocClin$sensitivities)
colnames(clinDF) <- c('Study', 'X', 'Y')

enetDF <- data.frame(rep('ElasticNet', length(rocEnet$specificities)),
                     1 - rocEnet$specificities, rocEnet$sensitivities)
colnames(enetDF) <- c('Study', 'X', 'Y')

pcrDF <- data.frame(rep('PrinCompReg', length(rocPcr$specificities)),
                    1 - rocPcr$specificities, rocPcr$sensitivities)
colnames(pcrDF) <- c('Study', 'X', 'Y')

plsDF <- data.frame(rep('PartLeastSq', length(rocPls$specificities)),
                    1 - rocPls$specificities, rocPls$sensitivities)
colnames(plsDF) <- c('Study', 'X', 'Y')

rfDF <- data.frame(rep('RandomForest', length(rocRF$specificities)),
                    1 - rocRF$specificities, rocRF$sensitivities)
colnames(rfDF) <- c('Study', 'X', 'Y')

fullDF <- rbind(clinDF, enetDF, pcrDF, plsDF, rfDF)

# Plot them
clinRocPlot <- ggplot(fullDF, aes(x = X, y = Y, group = Study)) + 
  geom_path(aes(colour = Study)) +
  geom_abline(slope = 1, colour = 'black') +
  opts(title = 'Receiver Operating Characteristic Curves\n') +
  xlab('1 - Specificity') +
  ylab('Sensitivity')
