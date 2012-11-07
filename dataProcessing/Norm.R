# Charles Ferté
# 09/07/2012
# Sage Bionetworks

###################################################################################
# Norm.R
# perform various unsupervised and supervised normalization methods on the different data
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

synapseLogin()

###################################################################################
# select dataset among Dir,Zhu, Hou, Lusc 
###################################################################################
dataset <- "Lusc"

###################################################################################
# load the data from Synapse
###################################################################################
Zhu <- loadEntity('syn1439020') #CEL files from Zhu et al
Dir <- loadEntity('syn1422422') #CEL Files from the Directors Challenge
Hou <- loadEntity('syn1422295') #CEL files from Hou et al
Lusc <- loadEntity('syn1426948') #CEL files from TCGA

###################################################################################
# now that the files are downloaded
# find them in the cache directory
###################################################################################
Dir_CEL <- Dir$cacheDir
Zhu_CEL <- Zhu$cacheDir
Hou_CEL <- Hou$cacheDir
Lusc_CEL <- Lusc$cacheDir

###################################################################################
# and now read the CEL files as affy objects
###################################################################################
datapath1  <-get(paste(dataset,"_CEL",sep=""))
setwd(datapath1)
rawdata <- ReadAffy(filenames=list.celfiles(datapath1))

##################################################################################
# make the sample Names coherent with the samples clinically curated
##################################################################################
# load the clinical data data from synapse
luscClin <- loadEntity('syn1438233')
luscClin <- luscClin$objects$LuscClinF

tmp <- intersect(sampleNames(rawdata), rownames(luscClin))
rawdata <- rawdata[,tmp]

###################################################################################
# perform rma normalization
###################################################################################
tmp <- paste(dataset,"_rma",sep="") #put rma at the end of the name to distinguish the results
assign(tmp,rma(rawdata))

###################################################################################
# perform gcrma normalization
###################################################################################
tmp <- paste(dataset,"_gcrma",sep="") #put gcrm at the end of the filenames
assign(tmp,gcrma(rawdata))

###################################################################################
# perform dCHIP normalization
###################################################################################
tmp <- paste(dataset,"_dCHIP",sep="") #put dCHIP at the end of the names
assign(tmp,expresso(rawdata, normalize.method = "invariantset", bg.correct = FALSE, pmcorrect.method = "pmonly",summary.method = "liwong"))

###################################################################################
# perform MAS5 normalization
###################################################################################
tmp <- paste(dataset,"_MAS5",sep="") #put MAS5 at the end of the names
assign(tmp,mas5(rawdata))

###################################################################################
# perform fRMA and then Barcode normalization
###################################################################################
tmp <- paste(dataset,"_frma",sep="") #put frma at the end of the names
assign(tmp,frma(rawdata, summarize = "random_effect"))
tmp1 <- paste(dataset,"_barcode",sep="")  #put barcode at the end of the names 
assign(tmp1,barcode(get(tmp)))

###################################################################################
# perform supervised normalization using snm
###################################################################################

# we know that GENDER and P_Stage are biological & study variables of interest 
bio.var <- model.matrix(~ luscClin$GENDER + luscClin$P_Stage)

# SCANBATCH is a variable concatenating the study name and the probable batch (grouped according to the cel files date of production)
adj.var <- model.matrix(~ luscClin$SITE )
myobject <- log2(pm(rawdata))
snm.fit <- snm(myobject, 
               bio.var=bio.var, 
               adj.var=adj.var, 
               rm.adj=TRUE)
new.expr <- snm.fit$norm.dat
pm(rawdata) <- 2^new.expr
myNormSummarized <- rma(rawdata, background=F, normalize=F)
dim(myNormSummarized)
tmp <- paste(dataset,"_snm",sep="") #put SNM at the end of the names
assign(tmp,myNormSummarized)

###################################################################################
# save the data in Synapse
#
# Saving steps for initial entities back to Synapse.  
# The ids will be accessible in the remainder of the code.
# all of this code is commented out to ensure that there is no change to the entities
# that already exist.
###################################################################################
#require(synapseClient)
#synapseLogin()
#
#dir_rma <- Data(list(name = "dir rma", parentId = 'syn87682'))
#dir_rma <- createEntity(dir_rma)

# add object into the data entity
#dir_rma <- addObject(dir_rma,Dir_rma)

# push the raw data into this entity
#dir_rma <- storeEntity(entity=dir_rma)

#############################################################################################
#dir_gcrma <- Data(list(name = "dir gcrma", parentId = 'syn87682'))
#dir_gcrma <- createEntity(dir_gcrma)

# add object into the data entity
#dir_gcrma <- addObject(dir_gcrma,Dir_gcrma)

# push the raw data into this entity
#dir_gcrma <- storeEntity(entity=dir_gcrma)
#############################################################################################


#dir_MAS5 <- Data(list(name = "dir MAS5", parentId = 'syn87682'))
#dir_MAS5 <- createEntity(dir_MAS5)

# add object into the data entity
#dir_MAS5 <- addObject(dir_MAS5,Dir_MAS5)

# push the raw data into this entity
#dir_MAS5 <- storeEntity(entity=dir_MAS5)
#############################################################################################
#dir_dCHIP <- Data(list(name = "dir dCHIP", parentId = 'syn87682'))
#dir_dCHIP <- createEntity(dir_dCHIP)

# add object into the data entity
#dir_dCHIP <- addObject(dir_dCHIP,Dir_dCHIP)

# push the raw data into this entity
#dir_dCHIP <- storeEntity(entity=dir_dCHIP)
#############################################################################################

#dir_frma <- Data(list(name = "dir frma", parentId = 'syn87682'))
#dir_frma <- createEntity(dir_frma)

# add object into the data entity
#dir_frma <- addObject(dir_frma,Dir_frma)

# push the raw data into this entity
#dir_frma <- storeEntity(entity=dir_frma)

#############################################################################################
# hou_snm <- Data(list(name = "hou_snm", parentId = 'syn87682'))
# hou_snm <- createEntity(hou_snm)
# 
# # add object into the data entity
# hou_snm <- addObject(hou_snm,Hou_snm)
# 
# # push the raw data into this entity
# hou_snm <- storeEntity(entity=hou_snm)


