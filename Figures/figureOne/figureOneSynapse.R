## figureOneSynapse.R

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washington
## erich.huang@sagebase.org

require(synapseClient)
## synapseLogin()

## FIGURE ONE PANEL ONE
figOnePanelOneEnt <- Data(list(name = 'figureOnePanelOneBinary',
                               parentId = 'syn1446518'))
figOnePanelOneEnt <- createEntity(figOnePanelOneEnt)
figOnePanelOneEnt <- addObject(figOnePanelOneEnt, rawPcPlot)
figOnePanelOneEnt <- storeEntity(figOnePanelOneEnt)