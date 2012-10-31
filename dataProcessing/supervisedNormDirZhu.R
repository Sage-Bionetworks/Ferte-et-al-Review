# Charles Ferté,MD MSc
# 10/11/2012
# Sage Bionetworks
#####
# code for supervised normalization of the Director's dataset (TS) and the Zhu dataset (VS)
#####

# load the packages that are neeeded
require(synapseClient)
require(affy)
require(caret)

synapseLogin()

######################################################################################################################################
# 1. load the datasets zhu and dir - preprocessing step to make them "comparable"
######################################################################################################################################

# load the rma normalized data from synapse
zhuExpr <- loadEntity('syn1436971')
zhuExpr <- zhuExpr$objects$Zhu_rma
dirExpr <- loadEntity('syn1440819')
dirExpr <- dirExpr$objects$Dir_rma

# load the clin data data from synapse
zhuClin <- loadEntity('syn1438225')
zhuClin <- zhuClin$objects$ZhuClinF
dirClin <- loadEntity('syn1438222')
dirClin <- dirClin$objects$DirClinF


# set P_Stage to be a factor with same levels across dir and zhu
dirClin$P_Stage <- factor(dirClin$P_Stage)
zhuClin$P_Stage <- factor(zhuClin$P_Stage, levels=levels(dirClin$P_Stage))


# Probability of 3 year overall survival is our dependant variable (vectors of response): os3yr
# make data coherent for expression and clinical data - only patients evaluable for os3yr
zhuClin$os3yr <-  ifelse(zhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36,1,ifelse(zhuClin$VITAL_STATUS==0,NA,0))
idzhu <- rownames(zhuClin)[ !is.na(zhuClin$os3yr) ]
zhuClin <- zhuClin[idzhu, ]
zhuExpr <- zhuExpr[, idzhu]
zhuExpr <- exprs(zhuExpr)

dirClin$os3yr <-  ifelse(dirClin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36,1,ifelse(dirClin$VITAL_STATUS==0,NA,0))
iddir <- rownames(dirClin)[ !is.na(dirClin$os3yr) ]
dirClin <- dirClin[iddir, ]
dirExpr <- dirExpr[, iddir]
dirExpr <- exprs(dirExpr)

# get rid of the control probesets
controlProbes <- grep("AFFX",rownames(zhuExpr))
zhuExpr <- zhuExpr[-controlProbes, ]
dirExpr <- dirExpr[rownames(zhuExpr), ]

# quantile normalize zhu to have the same disribution than dir (make them comparable)
zhuExpr <- normalize2Reference(zhuExpr, rowMeans(dirExpr))

# focus on the most variant probes
probVarZhu  <- apply(zhuExpr, 1, var)
probVarDir <- apply(dirExpr, 1, var)
meanProbVar <- apply(cbind(probVarZhu, probVarDir), 1, mean)
top20pct <- which(meanProbVar > quantile(meanProbVar, probs=.8))
dirExpr <- dirExpr[top20pct, ]
zhuExpr <- zhuExpr[top20pct, ]

# make the two datasets have the same mean and variance by scaling the validation data (zhuExpr)
# Justin Guinney's function to rescale the validation data to get the same mean/var than the training set
normalize_to_X <- function(mean.x, sd.x, Y){
  m.y <- rowMeans(Y)
  sd.y <- apply(Y, 1, sd)
  Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x 
  Y.adj[sd.y == 0] <- mean.x[sd.y==0]
  Y.adj
}

zhuExpr <- normalize_to_X(rowMeans(dirExpr), apply(dirExpr, 1, sd), zhuExpr)


####################################################################################################################################
# supervised normalization using snm adjusting for the inter-study batch, the intra batch, the histology, the stage and the gender
####################################################################################################################################
require(snm)
require(affy)

# load the raw data from dir and zhu
zhuraw <- loadEntity('syn1439020')
dirraw <- loadEntity('syn1422422')
filepath <- c(zhuraw$cacheDir,dirraw$cacheDir)
filenames <- list.celfiles(path=filepath,full.names=TRUE)
expr <- ReadAffy(filenames=filenames)
expr <- pm(expr)

# load the clinical data of dir and zhu
zhuclin <- loadEntity('syn1438225')
zhuclin <- zhuclin$objects$ZhuClinF
dirclin <- loadEntity('syn1438222')
dirclin <- dirclin$objects$DirClinF
allclin <- rbind(zhuclin,dirclin)

## compute the new data normalized 
# corrected for the SCANBATCH (batches determined by C.F according to the CEL file date of production and by the sudy status dir or zhu) 
bio.var <- model.matrix(~ allclin$GENDER + allclin$P_Stage)
adj.var <- model.matrix(~ allclin$SCANBATCH )
snm.fit <- snm(expr, 
               bio.var=bio.var, 
               adj.var=adj.var, 
               rm.adj=TRUE)

new.expr <- snm.fit$norm.dat



####################################################################################################################################


## PUSH INTERMEDIATE OBJECT TO SYNAPSE
supNormEnt <- Data(name="supervisedNormDirZhu", parentId="syn87682")
supNormEnt <- addObject(supNormEnt, zhuExpr)
supNormEnt <- addObject(supNormEnt, dirExpr)
supNormEnt <- addObject(supNormEnt, zhuClin)
supNormEnt <- addObject(supNormEnt, dirClin)
supNormEnt <- storeEntity(supNormEnt)
