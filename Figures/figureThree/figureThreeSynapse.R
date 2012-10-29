## figureThreeSynapse.R

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washington
## erich.huang@sagebase.org

## Save the figure binaries to Synapse
kmPlotList <- list('clinical' = clinKmPlot,
                   'elasticnet' = enetKmPlot,
                   'partleastsq' = plsKmPlot,
                   'princompreg' = pcrKmPlot,
                   'randomforest' = rfKmPlot)

kmPlotEnt <- Data(list(name = 'kmPlotListObject',
                       parentId = 'syn1449088'))
kmPlotEnt <- createEntity(kmPlotEnt)
kmPlotEnt <- addObject(kmPlotEnt, kmPlotList)
kmPlotEnt <- storeEntity(kmPlotEnt)

