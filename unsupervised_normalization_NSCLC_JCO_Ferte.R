# Charles Ferté
# 09/07/2012
#Sage Bionetworks

###################################################################################
# perform various unsupervised normalization methods on the different data
###################################################################################

###################################################################################
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
lol <- c(list.files(Dir_CEL,full.names=TRUE),list.files(Zhu_CEL,full.names=TRUE))
rawdata <- ReadAffy(filenames=as.character(lol))


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
HGU133a_rma <- Data(list(name = "HGU133a rma", parentId = 'syn87682'))
HGU133a_rma <- createEntity(HGU133a_rma)

# add object into the data entity
HGU133a_rma <- addObject(HGU133a_rma,HGU133A_rma)

# push the raw data into this entity
HGU133a_rma <- storeEntity(entity=HGU133a_rma)

#############################################################################################
HGU133a_gcrma <- Data(list(name = "HGU133a gcrma", parentId = 'syn87682'))
HGU133a_gcrma <- createEntity(HGU133a_gcrma)

# add object into the data entity
HGU133a_gcrma <- addObject(HGU133a_gcrma,HGU133A_gcrma)

# push the raw data into this entity
HGU133a_gcrma <- storeEntity(entity=HGU133a_gcrma)
#############################################################################################


HGU133a_MAS5 <- Data(list(name = "HGU133a MAS5", parentId = 'syn87682'))
HGU133a_MAS5 <- createEntity(HGU133a_MAS5)

# add object into the data entity
HGU133a_MAS5 <- addObject(HGU133a_MAS5,HGU133A_MAS5)

# push the raw data into this entity
HGU133a_MAS5 <- storeEntity(entity=HGU133a_MAS5)
#############################################################################################
HGU133a_dCHIP <- Data(list(name = "HGU133a dCHIP", parentId = 'syn87682'))
HGU133a_dCHIP <- createEntity(HGU133a_dCHIP)

# add object into the data entity
HGU133a_dCHIP <- addObject(HGU133a_dCHIP,HGU133A_dCHIP)

# push the raw data into this entity
HGU133a_dCHIP <- storeEntity(entity=HGU133a_dCHIP)
#############################################################################################

HGU133a_frma <- Data(list(name = "HGU133a_frma", parentId = 'syn87682'))
HGU133a_frma <- createEntity(HGU133a_frma)

# add object into the data entity
HGU133a_frma <- addObject(HGU133a_frma,HGU133A_frma)

# push the raw data into this entity
HGU133a_frma <- storeEntity(entity=HGU133a_frma)
#############################################################################################

HGU133a_barcode <- Data(list(name = "HGU133a barcode", parentId = 'syn87682'))
HGU133a_barcode <- createEntity(HGU133a_barcode)

# add object into the data entity
HGU133a_barcode <- addObject(HGU133a_barcode,HGU133A_barcode)

# push the raw data into this entity
HGU133a_barcode <- storeEntity(entity=HGU133a_barcode)
#
