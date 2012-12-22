# Charles Ferté,MD MSc
# 10/11/2012
# Sage Bionetworks

################################################################################################
# scalingDirZhu.R
# code for supervised normalization of the Director's dataset (TS) and the Zhu dataset (VS)
################################################################################################

# load the packages that are neeeded
require(synapseClient)
require(affy)
require(caret)

synapseLogin()

######################################################################################################################################
# 1. load the datasets zhu and dir - preprocessing step to make them "comparable"
######################################################################################################################################

# load the rma normalized data from synapse
zhuExpr <- loadEntity('syn1436971') #RMA normalized Zhu data
zhuExpr <- zhuExpr$objects$Zhu_rma
dirExpr <- loadEntity('syn1440819') #RMA normalized data from director's challenge
dirExpr <- dirExpr$objects$Dir_rma

# load the clinical data data from synapse
zhuClin <- loadEntity('syn1571258')
zhuClin <- zhuClin$objects$zhuClin
dirClin <- loadEntity('syn1571256')
dirClin <- dirClin$objects$dirClin

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

# describe the latent structure
vec <- c(rep(1,times=59),rep(2,times=274))
s <- svd(cbind(zhuExpr,dirExpr))
plot(s$v[,1],s$v[,2],col=c("royalblue","red")[vec],pch=20)

# focus on the most variant probes
probVarDir <- apply(dirExpr, 1, var)
top20pct <- which(probVarDir > quantile(probVarDir, probs=.8))
dirExpr <- dirExpr[top20pct, ]
zhuExpr <- zhuExpr[top20pct, ]

# describe the latent structure
s <- svd(cbind(zhuExpr,dirExpr))
plot(s$v[,1],s$v[,2],col=c("royalblue","red")[vec],pch=20)

# make the two datasets have the same mean and variance by scaling the validation data (zhuExpr) using scale function
zhuExpr1 <- t(scale(t(zhuExpr),center=TRUE,scale=TRUE))
dirExpr1 <- t(scale(t(dirExpr),center=TRUE,scale=TRUE))



# Justin Guinney's function to rescale the validation data to get the same mean/var than the training set
normalize_to_X <- function(mean.x, sd.x, Y){
 m.y <- rowMeans(Y)
 sd.y <- apply(Y, 1, sd)
 Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x 
 Y.adj[sd.y == 0] <- mean.x[sd.y==0]
 Y.adj
}

dirExpr1 <- dirExpr
zhuExpr1 <- normalize_to_X(rowMeans(dirExpr1), apply(dirExpr1, 1, sd), zhuExpr)

# describe the latent structure
s <- svd(cbind(zhuExpr1,dirExpr1))
plot(s$v[,1],s$v[,2],col=c("royalblue","red")[vec],pch=20)

dirExpr <- dirExpr1
zhuExpr <- zhuExpr1

###################################################################################
# save the data in Synapse
#
# Saving steps for initial entities back to Synapse.  
# The ids will be accessible in the remainder of the code.
# all of this code is commented out to ensure that there is no change to the entities
# that already exist.
###################################################################################

# supNormEnt <- Data(name="supervisedNormDirZhu", parentId="syn87682")
# supNormEnt <- addObject(supNormEnt, zhuExpr)
# supNormEnt <- addObject(supNormEnt, dirExpr)
# supNormEnt <- addObject(supNormEnt, zhuClin)
# supNormEnt <- addObject(supNormEnt, dirClin)
# supNormEnt <- storeEntity(supNormEnt)


# scaledfRMAEnt <- Data(name="scaledfRMADirZhu", parentId="syn87682")
# scaledfRMAEnt <- addObject(scaledfRMAEnt, zhuExpr)
# scaledfRMAEnt <- addObject(scaledfRMAEnt, dirExpr)
# scaledfRMAEnt <- addObject(scaledfRMAEnt, zhuClin)
# scaledfRMAEnt <- addObject(scaledfRMAEnt, dirClin)
# scaledfRMAEnt <- storeEntity(scaledfRMAEnt)

# scaledBarcodeEnt <- Data(name="scaledBarcodeDirZhu", parentId="syn87682")
# scaledBarcodeEnt <- addObject(scaledBarcodeEnt, zhuExpr)
# scaledBarcodeEnt <- addObject(scaledBarcodeEnt, dirExpr)
# scaledBarcodeEnt <- addObject(scaledBarcodeEnt, zhuClin)
# scaledBarcodeEnt <- addObject(scaledBarcodeEnt, dirClin)
# scaledBarcodeEnt <- storeEntity(scaledBarcodeEnt)
