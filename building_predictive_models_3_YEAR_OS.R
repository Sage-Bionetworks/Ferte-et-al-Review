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
zhu_clin$y_zhu <- ifelse(zhu_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36,1,0)
dir_clin$y_dir <- ifelse(dir_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36,1,0)

# first build the model based on clin variables of interest only (pStage, gender, age,smoking) (linear regression)
fit <- glm(dir_clin$y_dir~dir_clin$P_Stage + dir_clin$GENDER + dir_clin$Age,data=as.data.frame(t(dir_clin)),family="binomial")
summary(fit)
confint.default(fit)
boxplot(fit$fitted.values~dir_clin$y_dir,ylab="3-year OS prediction (%)",xlab="3-year OS")
stripchart(fit$fitted.values~dir_clin$y_dir,pch=20,col="royalblue",vertical=TRUE,add=TRUE,cex=.6)

# predict in zhu
yhat_zhu <- predict.glm(fit,newdata=as.data.frame(zhu_clin),type="response",na.action = na.omit)

# second build the model based on linear regression of clin + ge features


# third build the model based on molecular features using elasticnet
cv.fit <- cv.glmnet(t(exprs(dir)),y=y_dir,nfolds=10,alpha=.1,family="binomial")
plot(cv.fit)
fit <- glmnet(family="binomial",alpha=.1,lambda=cv.fit$lambda.min,y=y_dir,x=t(exprs(dir)))

yhat <- predict()
# thrird build the model based on clin + molecular features using RandomForest

# adress the results of in the VS using ROC curves

# adress the same results of AUC in boxplots with their IC95

