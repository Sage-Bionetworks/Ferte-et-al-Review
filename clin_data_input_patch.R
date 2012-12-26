# charles ferté
# dec 22 2012


# correct the clin files of dir zhu lusc and hou
# problem was the "THREE_YEAR_OS"...which was wrong 


# Load the data 
houClin <- loadEntity('syn1438227')
houClin <- houClin$objects$HouClinF

# Load the data 
luscClin <- loadEntity('syn1438233')
luscClin <- luscClin$objects$LuscClinF

# correct the error by removing this variable
dirClin$THREE_YEAR_OS <- NULL
zhuClin$THREE_YEAR_OS <- NULL
houClin$THREE_YEAR_OS <- NULL
luscClin$THREE_YEAR_OS <- NULL

# save the data under Synapse entity syn1488297
LuscClin <- Data(list(name = "LuscClin", parentId = 'syn1488297'))
LuscClin <- createEntity(LuscClin)
# add object into the data entity
LuscClin <- addObject(LuscClin,luscClin)
# push the raw data into this entity
LuscClin <- storeEntity(entity=LuscClin)





