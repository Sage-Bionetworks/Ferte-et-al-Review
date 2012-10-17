# Charles Ferté
# 09/07/2012
# Sage Bionetworks

###################################################################################
# perform various unsupervised normalization methods on the different data
###################################################################################

##################################################################################
# call the necessary libraries
###################################################################################
library(affy)
library(gcrma)
library(frma)
library(hgu133plus2frmavecs)
library(hgu133afrmavecs)
library(hthgu133afrmavecs)
library(snm)
require(synapseClient)
synapseLogin(username="charles.ferte@sagebase.org",password="charles")

###################################################################################
# select dataset among Dir,Zhu, Hou, Lusc or HGU133A (aka combined Dir and Zhu HGU133A arrays)
###################################################################################
dataset <- "Dir"

###################################################################################
# load the data from Synapse
###################################################################################
Zhu <- loadEntity('syn1421817')
Dir <- loadEntity('syn1422422')
Hou <- loadEntity('syn1422295')
Lusc <- loadEntity('syn1426948')

###################################################################################
# where are the CEL files
###################################################################################
Dir_CEL <- Dir$cacheDir
Zhu_CEL <- Zhu$cacheDir
Hou_CEL <- Hou$cacheDir
Lusc_CEL <- Lusc$cacheDir

###################################################################################
# read the CEL files as affy objects
###################################################################################
datapath1  <-get(paste(dataset,"_CEL",sep=""))
setwd(datapath1)
rawdata <- ReadAffy()

# read the HGU133A
#lol <- c(list.files(Dir_CEL,full.names=TRUE),list.files(Zhu_CEL,full.names=TRUE))
#rawdata <- ReadAffy(filenames=as.character(lol))


###################################################################################
# perform rma normalization
###################################################################################
tmp <- paste(dataset,"_rma",sep="")
assign(tmp,rma(rawdata))

###################################################################################
# perform gcrma normalization
###################################################################################
tmp <- paste(dataset,"_gcrma",sep="")
assign(tmp,gcrma(rawdata))

###################################################################################
# perform dCHIP normalization
###################################################################################
tmp <- paste(dataset,"_dCHIP",sep="")
assign(tmp,expresso(rawdata, normalize.method = "invariantset", bg.correct = FALSE, pmcorrect.method = "pmonly",summary.method = "liwong"))

###################################################################################
# perform MAS5 normalization
###################################################################################
tmp <- paste(dataset,"_MAS5",sep="")
assign(tmp,mas5(rawdata))

###################################################################################
# perform fRMA and then Barcode normalization
###################################################################################
tmp <- paste(dataset,"_frma",sep="")
assign(tmp,frma(rawdata, summarize = "random_effect"))
tmp1 <- paste(dataset,"_barcode",sep="")
assign(tmp1,barcode(get(tmp)))

###################################################################################
# save the data in Synapse
###################################################################################
require(synapseClient)
synapseLogin(username="charles.ferte@sagebase.org",password="charles")
#
dir_rma <- Data(list(name = "dir rma", parentId = 'syn87682'))
dir_rma <- createEntity(dir_rma)

# add object into the data entity
dir_rma <- addObject(dir_rma,Dir_rma)

# push the raw data into this entity
dir_rma <- storeEntity(entity=dir_rma)

#############################################################################################
dir_gcrma <- Data(list(name = "dir gcrma", parentId = 'syn87682'))
dir_gcrma <- createEntity(dir_gcrma)

# add object into the data entity
dir_gcrma <- addObject(dir_gcrma,Dir_gcrma)

# push the raw data into this entity
dir_gcrma <- storeEntity(entity=dir_gcrma)
#############################################################################################


dir_MAS5 <- Data(list(name = "dir MAS5", parentId = 'syn87682'))
dir_MAS5 <- createEntity(dir_MAS5)

# add object into the data entity
dir_MAS5 <- addObject(dir_MAS5,Dir_MAS5)

# push the raw data into this entity
dir_MAS5 <- storeEntity(entity=dir_MAS5)
#############################################################################################
dir_dCHIP <- Data(list(name = "dir dCHIP", parentId = 'syn87682'))
dir_dCHIP <- createEntity(dir_dCHIP)

# add object into the data entity
dir_dCHIP <- addObject(dir_dCHIP,Dir_dCHIP)

# push the raw data into this entity
dir_dCHIP <- storeEntity(entity=dir_dCHIP)
#############################################################################################

dir_frma <- Data(list(name = "dir frma", parentId = 'syn87682'))
dir_frma <- createEntity(dir_frma)

# add object into the data entity
dir_frma <- addObject(dir_frma,Dir_frma)

# push the raw data into this entity
dir_frma <- storeEntity(entity=dir_frma)
#############################################################################################

dir_barcode <- Data(list(name = "dir barcode", parentId = 'syn87682'))
dir_barcode <- createEntity(dir_barcode)

# add object into the data entity
dir_barcode <- addObject(dir_barcode,Dir_barcode)

# push the raw data into this entity
dir_barcode <- storeEntity(entity=dir_barcode)
#
