## figureOneSynapse.R

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washington
## erich.huang@sagebase.org

require(synapseClient)
## synapseLogin()

## FIGURE HELPER FUNCTIONS
figFunctionsEnt <- Code(list(name = 'figureFunctions.R',
                        parentId = 'syn87682'))
figFunctionsEnt <- createEntity(figFunctionsEnt)
figFunctionsEnt <- addFile(figFunctionsEnt, 'figureFunctions.R')
figFunctionsEnt <- storeEntity(figFunctionsEnt)

## FIGURE ONE PANEL ONE
figOnePanelOneEnt <- Data(list(name = 'figureOnePanelOneBinary',
                               parentId = 'syn1446518'))
figOnePanelOneEnt <- createEntity(figOnePanelOneEnt)
figOnePanelOneEnt <- addObject(figOnePanelOneEnt, rawPcPlot)
figOnePanelOneEnt <- storeEntity(figOnePanelOneEnt)

## FIGURE ONE PANEL TWO
figOnePanelTwoEnt <- Data(list(name = 'figureOnePanelTwoBinary',
                               parentId = 'syn1446518'))
figOnePanelTwoEnt <- createEntity(figOnePanelTwoEnt)
figOnePanelTwoEnt <- addObject(figOnePanelTwoEnt, mas5PcPlot)
figOnePanelTwoEnt <- storeEntity(figOnePanelTwoEnt)

## FIGURE ONE PANEL THREE
figOnePanelThreeEnt <- Data(list(name = 'figureOnePanelThreeBinary',
                                 parentId = 'syn1446518'))
figOnePanelThreeEnt <- createEntity(figOnePanelThreeEnt)
figOnePanelThreeEnt <- addObject(figOnePanelThreeEnt, rmaPcPlot)
figOnePanelThreeEnt <- storeEntity(figOnePanelThreeEnt)

## FIGURE ONE PANEL FOUR
figOnePanelFourEnt <- Data(list(name = 'figureOnePanelFourBinary',
                                parentId = 'syn1446518'))
figOnePanelFourEnt <- createEntity(figOnePanelFourEnt)
figOnePanelFourEnt <- addObject(figOnePanelFourEnt, gcrmaPcPlot)
figOnePanelFourEnt <- storeEntity(figOnePanelFourEnt)

## FIGURE ONE PANEL FIVE
figOnePanelFiveEnt <- Data(list(name = 'figureOnePanelFiveBinary',
                                parentId = 'syn1446518'))
figOnePanelFiveEnt <- createEntity(figOnePanelFiveEnt)
figOnePanelFiveEnt <- addObject(figOnePanelFiveEnt, dchipPcPlot)
figOnePanelFiveEnt <- storeEntity(figOnePanelFiveEnt)


