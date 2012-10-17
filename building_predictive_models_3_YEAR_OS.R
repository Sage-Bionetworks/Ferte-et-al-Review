# Charles Ferté
# 10/11/2012
#Sage Bionetworks

# code to develop the model in the Director's dataset (TS) and validate it in the Zhu dataset (VS)
# we want to predict the 3 year propbability of OS

# load the packages that are neeeded
require(glmnet)
require(randomForest)
require(synapseClient)
require(affy)
synapseLogin(username="charles.ferte@sagebase.org",password="charles")

# load the rma normalized data from synapse
zhu <- loadEntity('syn1436971')
zhu <- zhu$objects$Zhu_rma
dir <- loadEntity('syn1440819')
dir <- dir$objects$Dir_rma

# load the clin data data from synapse
zhu_clin <- loadEntity('syn1438225')
zhu_clin <- zhu_clin$objects$ZhuClinF
dir_clin <- loadEntity('syn1438222')
dir_clin <- dir_clin$objects$DirClinF

# make data coherent for dir since there are two cel files that are not clinically annotated
tmp <- intersect(rownames(dir_clin),sampleNames(dir))
dir_clin <- dir_clin[tmp,]
dir <- dir[,tmp]
rm(tmp)

# create vectors of response( y_dir and y_zhu) of 3 year overall survival (1 alive, 0 deceased)
y_zhu <- ifelse(zhu_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36,1,0)
names(y_zhu) <- rownames(zhu_clin)

y_dir <- ifelse(dir_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36,1,0)
names(y_dir) <- rownames(dir_clin)

# first build the model based on clin variables only (linear regression)

# second build the model based on linear regression of clin + ge features

# thrird build the model based on clin + molecular features using elasticnet

# thrird build the model based on clin + molecular features using RandomForest

# adress the results of in the VS using ROC curves

# adress the same results of AUC in boxplots with their IC95

