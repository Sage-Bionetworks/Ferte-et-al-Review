# Charles Ferté
# 09/07/2012
# Sage Bionetworks

###################################################################################
# perform various unsupervised normalization methods
###################################################################################

##################################################################################
# call the required libraries
###################################################################################
require(affy)
require(gcrma)
require(frma)
require(snm)
require(synapseClient)
require(hgu133plus2frmavecs)
require(hgu133afrmavecs)
require(hthgu133afrmavecs)

synapseLogin()

###################################################################################
# load the data from Synapse
###################################################################################
Dir <- loadEntity('syn1422422')

# CEL file directory
Dir_CEL <- Dir$cacheDir

# read the CEL files as affy objects
dataPath  <- "Dir_CEL"
rawdata <- ReadAffy(filenames=list.celfiles(dataPath))

# perform rma normalization
Dir_rma <- rma(rawdata)

# perform gcrma normalization
Dir_gcrma <- gcrma(rawdata)

# perform dCHIP normalization
Dir_dCHIP <- expresso(rawdata, normalize.method = "invariantset", bg.correct = FALSE, pmcorrect.method = "pmonly", summary.method = "liwong")

# perform MAS5 normalization
Dir_MAS5 <- mas5(rawdata)

# perform fRMA and then Barcode normalization
Dir_frma <- frma(rawdata, summarize = "random_effect")
Dir_barcode <- barcode(Dir_frma)

###################################################################################
# save the data in Synapse
#
# Saving steps for initial entities back to Synapse.  
# The ids will be accessible in the remainder of the code.
###################################################################################

#############################################################################################
## rma
dir_rma <- Data(list(name = "dir rma", parentId = 'syn87682'))
dir_rma <- createEntity(dir_rma)

# add object into the data entity
dir_rma <- addObject(dir_rma,Dir_rma)

# push the raw data into this entity
dir_rma <- storeEntity(dir_rma)

#############################################################################################
## gcrma
dir_gcrma <- Data(list(name = "dir gcrma", parentId = 'syn87682'))
dir_gcrma <- createEntity(dir_gcrma)

# add object into the data entity
dir_gcrma <- addObject(dir_gcrma,Dir_gcrma)

# push the raw data into this entity
dir_gcrma <- storeEntity(dir_gcrma)

#############################################################################################
## MAS5
dir_MAS5 <- Data(list(name = "dir MAS5", parentId = 'syn87682'))
dir_MAS5 <- createEntity(dir_MAS5)

# add object into the data entity
dir_MAS5 <- addObject(dir_MAS5,Dir_MAS5)

# push the raw data into this entity
dir_MAS5 <- storeEntity(dir_MAS5)

#############################################################################################
## dCHIP
dir_dCHIP <- Data(list(name = "dir dCHIP", parentId = 'syn87682'))
dir_dCHIP <- createEntity(dir_dCHIP)

# add object into the data entity
dir_dCHIP <- addObject(dir_dCHIP,Dir_dCHIP)

# push the raw data into this entity
dir_dCHIP <- storeEntity(dir_dCHIP)

#############################################################################################
## frma
dir_frma <- Data(list(name = "dir frma", parentId = 'syn87682'))
dir_frma <- createEntity(dir_frma)

# add object into the data entity
dir_frma <- addObject(dir_frma,Dir_frma)

# push the raw data into this entity
dir_frma <- storeEntity(dir_frma)

#############################################################################################
## barcode
dir_barcode <- Data(list(name = "dir barcode", parentId = 'syn87682'))
dir_barcode <- createEntity(dir_barcode)

# add object into the data entity
dir_barcode <- addObject(dir_barcode,Dir_barcode)

# push the raw data into this entity
dir_barcode <- storeEntity(dir_barcode)

