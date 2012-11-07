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
# perform fRMA normalization
###################################################################################
tmp <- paste(dataset,"_frma",sep="") #put frma at the end of the names
assign(tmp,frma(rawdata, summarize = "random_effect"))

###################################################################################
# perform supervised normalization using snm (and rma summarization)
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

##################################################################################
# save the data in Synapse
# 
# Saving steps for initial entities back to Synapse.  
# The ids will be accessible in the remainder of the code.
# all of this code is commented out to ensure that there is no change to the entities
# that already exist.
##################################################################################
require(synapseClient)
synapseLogin()

lusc_rma <- Data(list(name = "lusc_rma", parentId = 'syn87682'))
lusc_rma <- createEntity(lusc_rma)

# add object into the data entity
lusc_rma <- addObject(lusc_rma,Lusc_rma)

# push the raw data into this entity
lusc_rma <- storeEntity(entity=lusc_rma)

############################################################################################

lusc_gcrma <- Data(list(name = "lusc_gcrma", parentId = 'syn87682'))
lusc_gcrma <- createEntity(lusc_gcrma)

# add object into the data entity
lusc_gcrma <- addObject(lusc_gcrma,Lusc_gcrma)

# push the raw data into this entity
lusc_gcrma <- storeEntity(entity=lusc_gcrma)

############################################################################################

lusc_MAS5 <- Data(list(name = "lusc_MAS5", parentId = 'syn87682'))
lusc_MAS5 <- createEntity(lusc_MAS5)

#add object into the data entity
lusc_MAS5 <- addObject(lusc_MAS5,Lusc_MAS5)

#push the raw data into this entity
lusc_MAS5 <- storeEntity(entity=lusc_MAS5)

############################################################################################

lusc_dCHIP <- Data(list(name = "lusc_dCHIP", parentId = 'syn87682'))
lusc_dCHIP <- createEntity(lusc_dCHIP)

# add object into the data entity
lusc_dCHIP <- addObject(lusc_dCHIP,Lusc_dCHIP)

# push the raw data into this entity
lusc_dCHIP <- storeEntity(entity=lusc_dCHIP)

############################################################################################

lusc_frma <- Data(list(name = "lusc_frma", parentId = 'syn87682'))
lusc_frma <- createEntity(lusc_frma)

#add object into the data entity
lusc_frma <- addObject(lusc_frma,Lusc_frma)

#push the raw data into this entity
lusc_frma <- storeEntity(entity=lusc_frma)

############################################################################################
lusc_snm <- Data(list(name = "lusc_snm", parentId = 'syn87682'))
lusc_snm <- createEntity(lusc_snm)

# add object into the data entity
lusc_snm <- addObject(lusc_snm,Lusc_snm)

# push the raw data into this entity
lusc_snm <- storeEntity(entity=lusc_snm)


