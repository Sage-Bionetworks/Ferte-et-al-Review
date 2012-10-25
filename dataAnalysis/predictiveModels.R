# Charles Ferté,MD MSc
# 10/11/2012
# Sage Bionetworks

# code to develop the model in the Director's dataset (TS) and validate it in the Zhu dataset (VS)
# the dependant variable (i.e. variable to be predicted) is the 3 year propbability of OS


# load the packages that are neeeded
require(glmnet)
require(randomForest)
require(synapseClient)
require(caret)
require(affy)
require(survival)
require(ROCR)
require(pls)
set.seed(12221981)
######################################################################################################################################
# 1. load the datasets zhu and dir - preprocessing step to make them "comparable"
######################################################################################################################################

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

# Probability of 3 year overall survival is our dependant variable (vectors of response): y_dir and y_zhu
# make data coherent for expression and clinical data
zhu_clin$y_zhu <-  ifelse(zhu_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36,1,ifelse(zhu_clin$VITAL_STATUS==0,NA,0))
idzhu <- rownames(zhu_clin)[ !is.na(zhu_clin$y_zhu) ]
zhu_clin <- zhu_clin[idzhu, ]
zhu <- zhu[, idzhu]
zhu <- exprs(zhu)

dir_clin$y_dir <-  ifelse(dir_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36,1,ifelse(dir_clin$VITAL_STATUS==0,NA,0))
iddir <- rownames(dir_clin)[ !is.na(dir_clin$y_dir) ]
dir_clin <- dir_clin[iddir, ]
dir <- dir[, iddir]
dir <- exprs(dir)

# get rid of the control probesets
controlProbes <- grep("AFFX",rownames(zhu))
zhu <- zhu[-controlProbes, ]
dir <- dir[rownames(zhu), ]

# quantile normalize zhu to have the same disribution than dir (make them comparable)
zhu <- normalize2Reference(zhu, rowMeans(dir))

# focus on the most variant probes
prob.var1  <- apply(zhu, 1, var)
prob.var2 <- apply(dir, 1, var)
mean.prob.var <- apply(cbind(prob.var1, prob.var2), 1, mean)
tmp <- which(mean.prob.var > quantile(mean.prob.var, probs=.8))
dir <- dir[tmp, ]
zhu <- zhu[tmp, ]
rm(tmp, mean.prob.var, prob.var1, prob.var2)

# make the two datasets have the same mean and variance
# Justin Guinney's function to rescale the validation data to get the same mean/var than the training set
normalize_to_X <- function(mean.x, sd.x, Y){
  m.y <- rowMeans(Y)
  sd.y <- apply(Y, 1, sd)
  Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x 
  Y.adj[sd.y == 0] <- mean.x[sd.y==0]
  Y.adj
}

zhu <- normalize_to_X(rowMeans(dir),apply(dir,1,sd),zhu)

######################################################################################################################################
# 2. build a model based on clin variables of interest only (pStage, gender, age) (logisitic regression)
######################################################################################################################################

dir_clin$P_Stage <- factor(dir_clin$P_Stage)
zhu_clin$P_Stage <- factor(zhu_clin$P_Stage, levels=levels(dir_clin$P_Stage))

fitClin <- glm(y_dir ~ P_Stage, data = dir_clin, family = "binomial")
#summary(fitClin)
yhatClin <- predict(fitClin, zhu_clin, type = "response")

#########################################################################################################################################
# 3. build the model based on molecular features using elasticnet 
#  not penalizing the clinical features in a more ridge setting (alpha=.1) 
#########################################################################################################################################

# let x and z be the datsets combining molecular features and P Stage
dirP <- t(model.matrix(~ -1 + dir_clin$P_Stage))
rownames(dirP) <- sub("dir_clin$", "", rownames(dirP), fixed=T)
rownames(dirP) <- sub(")", "_", rownames(dirP), fixed=t)
x <- rbind(dir, dirP)

zhuP <- t(model.matrix(~ -1 + zhu_clin$P_Stage))
rownames(zhuP) <- sub("zhu_clin$", "", rownames(zhuP), fixed=T)
rownames(zhuP) <- sub(")", "_", rownames(zhuP), fixed=t)
z <- rbind(zhu, zhuP)

# let pen be a vector of penalty factors for each coefficient. We want to preserve the pathological stage in the model, so we assign it 0 as penalty factor.
pen <- c(rep(1, nrow(dir)), rep(0, nrow(dirP)))

# run the elastic net model
cv.fit <- cv.glmnet(x=t(x), y=factor(dir_clin$y_dir), nfolds=10, alpha=.1, family="binomial", penalty.factor=pen)
plot(cv.fit)
fitEnet <- glmnet(x=t(x), y=factor(dir_clin$y_dir), family="binomial", alpha=.1, lambda=cv.fit$lambda.min, penalty.factor=pen)
yhatEnet <- predict(fitEnet, t(z), type="response", s="lambda.min")

boxplot(yhatEnet ~ zhu_clin$y_zhu, ylab="3-year OS prediction (%)", xlab="3-year OS", main="elastic net - molecular + clinical features")
stripchart(yhatEnet ~ zhu_clin$y_zhu,pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

######################################################################################################################################
# 4. build the model based on clin + molecular features using RandomForest
######################################################################################################################################

# let's find the optimal value of mtry
a <- tuneRF(x=t(x), y=factor(dir_clin$y_dir), stepFactor=1.5, doBest=TRUE, trace=TRUE, ntreeTry=1000)

# run the model with this mtry value
fitRF <- randomForest(x=t(x), y=factor(dir_clin$y_dir), mtry=a$mtry, do.trace=10, ntree=1000, importance=TRUE)
yhatRF <- predict(fitRF, t(z), type="prob")[,2]

boxplot(yhatRF ~ zhu_clin$y_zhu, ylab="3-year OS prediction (%)", xlab="3-year OS", main= "Random Forest - molecular + clinical features")
stripchart(yhatRF ~ zhu_clin$y_zhu, pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

###################################################################################################################
# 5. principal component regression
##############################################################################################################

fitPcr <- pcr(dir_clin$y_dir ~ t(x), ncomp=10,validation = "CV", family="binomial")
yhatPcr <- predict(fitPcr, comps = 1:9,t(z), type="response")

boxplot(yhatPcr ~ zhu_clin$y_zhu, ylab="3-year OS prediction (%)", xlab="3-year OS", main= "Principal component regression - molecular + clinical features")
stripchart(yhatPcr ~ zhu_clin$y_zhu, pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

###################################################################################################################
# 6. partial leased square
##############################################################################################################
fitPls <- plsr(dir_clin$y_dir ~ t(x), ncomp=10,validation = "CV", family="binomial")
yhatPls <- predict(fitPls, comps = 1:9,t(z), type="response")

boxplot(yhatPls ~ zhu_clin$y_zhu, ylab="3-year OS prediction (%)", xlab="3-year OS", main= "Principal component regression - molecular + clinical features")
stripchart(yhatPls ~ zhu_clin$y_zhu, pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

######################################################################################################################################
# 7. plot a ROC curve to asses the performance of our models
######################################################################################################################################

Pred <- prediction(as.numeric(yhatClin), as.numeric(zhu_clin$y_zhu))
Perf <- performance(prediction.obj=Pred, "tpr", "fpr")
AUC <- performance(prediction.obj=Pred, "auc")
plot(Perf, col="royalblue", main="performance of the models 
predicting the probability of 3 years OS", lwd=2)
text(x=.25, y=.25, labels=paste("AUC logit clin. =", format(x=AUC@y.values, digits=2)), col="royalblue", adj=0, cex=.8)

Pred <- prediction(as.numeric(yhatEnet), as.numeric(zhu_clin$y_zhu))
Perf <- performance(prediction.obj=Pred, "tpr", "fpr")
AUC <- performance(prediction.obj=Pred, "auc")
plot(Perf, col="red", lwd=2, add=TRUE)
text(x=.25, y=.2, labels=paste("AUC Elastic Net=", format(x=AUC@y.values, digits=2)), col="red", adj=0, cex=.8)

Pred <- prediction(as.numeric(yhatRF), as.numeric(zhu_clin$y_zhu))
Perf <- performance(prediction.obj=Pred, "tpr", "fpr")
AUC <- performance(prediction.obj=Pred, "auc")
plot(Perf, col="orange", add=TRUE, lwd=2)
text(x=.25, y=.15, labels=paste("AUC Random Forest=", format(x=AUC@y.values, digits=2)), col="orange", adj=0, cex=.8)

Pred <- prediction(as.numeric(yhatPcr), as.numeric(zhu_clin$y_zhu))
Perf <- performance(prediction.obj=Pred, "tpr", "fpr")
AUC <- performance(prediction.obj=Pred, "auc")
plot(Perf, col="darkgreen", add=TRUE, lwd=2)
text(x=.25, y=.1, labels=paste("AUC Principal Comp. regression=", format(x=AUC@y.values, digits=2)), col="darkgreen", adj=0, cex=.8)

Pred <- prediction(as.numeric(yhatPls), as.numeric(zhu_clin$y_zhu))
Perf <- performance(prediction.obj=Pred, "tpr", "fpr")
AUC <- performance(prediction.obj=Pred, "auc")
plot(Perf, col="black", add=TRUE, lwd=2)
text(x=.25, y=.05, labels=paste("AUC Partial least square=", format(x=AUC@y.values, digits=2)), col="black", adj=0, cex=.8)


######################################################################################################################################
# 7. draw the kaplan meier curves based on the predictors (high and low risk groups based on the median)
######################################################################################################################################

# for each model, let's assign to "high-risk" the group of patients whose yhat > median(yhat) 
# and to "low-risk" the group of patients whose yhat<median(yhat). We set high risk =1, low risk =0.

riskClin <- ifelse(yhatClin >= median(yhatClin), 1, 0)
names(riskClin) <- names(yhatClin)
riskEnet <- ifelse(yhatEnet >= median(yhatEnet), 1, 0)
names(riskEnet) <- names(yhatEnet)
riskRF <- ifelse(yhatRF >= median(yhatRF), 1, 0)
riskPcr <- ifelse(yhatPcr >= median(yhatPcr), 1, 0)
riskPls <- ifelse(yhatPls >= median(yhatPls), 1, 0)

# for each model, draw the kaplan meir curves and compute the log rank test

par(mfrow=c(2,3))
plot(survfit(Surv(zhu_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhu_clin$VITAL_STATUS) ~ riskClin), main="logit model", xlab="months",ylab="probability of 3 years OS")
survdiff(Surv(zhu_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhu_clin$VITAL_STATUS) ~ riskClin, rho=0)

plot(survfit(Surv(zhu_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhu_clin$VITAL_STATUS) ~ riskEnet), main="elastic net model", xlab="months",ylab="probability of 3 years OS")
survdiff(Surv(zhu_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhu_clin$VITAL_STATUS) ~ riskEnet, rho=0)

plot(survfit(Surv(zhu_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhu_clin$VITAL_STATUS) ~ riskRF), main="random forest model", xlab="months",ylab="probability of 3 years OS")
survdiff(Surv(zhu_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhu_clin$VITAL_STATUS) ~ riskRF, rho=0)

plot(survfit(Surv(zhu_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhu_clin$VITAL_STATUS) ~ riskPcr), main="Prin. Comp. Regression model", xlab="months",ylab="probability of 3 years OS")
survdiff(Surv(zhu_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhu_clin$VITAL_STATUS) ~ riskPcr, rho=0)

plot(survfit(Surv(zhu_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhu_clin$VITAL_STATUS) ~ riskPls), main="Partial least square model", xlab="months",ylab="probability of 3 years OS")
survdiff(Surv(zhu_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhu_clin$VITAL_STATUS) ~ riskPls, rho=0)


######################################################################################################################################
# 7. draw the heatmaps based on the predictors (high and low risk groups based on the median)
######################################################################################################################################
require(gplots)
sort(riskEnet)
orderCol <- rownames(riskEnet)[sort(riskEnet,index.return=T)$ix]
heatmap.2(x=z[which(abs(fitEnet$beta)>0),orderCol],trace="none",ColSideColors=c("blue","red")[riskEnet[orderCol]+1],col=greenred(100),Colv=FALSE)
          
dim(x[which(abs(fitEnet$beta)>0),])