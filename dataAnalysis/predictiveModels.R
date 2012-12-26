# Charles Ferté,MD MSc
# 10/11/2012
# Sage Bionetworks


# code to develop the model in the Director's dataset (TS) and validate it in the Zhu dataset (VS)
# the dependant variable (i.e. variable to be predicted) is the 3 year propbability of OS

# load the packages that are neeeded
require(synapseClient)
require(glmnet)
require(randomForest)
require(caret)
require(affy)
require(survival)
require(pROC)
require(pls)
require(gplots)

set.seed(12221981)

synapseLogin()

######################################################################################################################################
# 1. load the preprocessed datasets zhu and dir
######################################################################################################################################

# load the supervised normalized data from synapse
supNormEnt <- loadEntity("syn1571600")
dirClin <- supNormEnt$objects$dirClin
dirExpr <- supNormEnt$objects$dirExpr
zhuClin <- supNormEnt$objects$zhuClin
zhuExpr <- supNormEnt$objects$zhuExpr


######################################################################################################################################
# 2. build a model of os3yr based on P_Stage only since this is the current clinical standard (logisitic regression)
######################################################################################################################################

fitClin <- glm(os3yr ~ P_Stage + GENDER, data = dirClin, family = "binomial")
summary(fitClin)

yhatClin <- predict(fitClin, zhuClin, type = "response")

#########################################################################################################################################
# 3. build the model based on molecular features using elasticnet 
#  not penalizing the clinical features in a more ridge setting (alpha=.1) 
#########################################################################################################################################

# let x and z be the datsets combining molecular features and P Stage
dirP <- t(model.matrix(~ -1 + dirClin$P_Stage))
rownames(dirP) <- sub("dirClin$", "", rownames(dirP), fixed=T)
rownames(dirP) <- sub(")", "_", rownames(dirP), fixed=t)
x <- rbind(dirExpr, dirP)

zhuP <- t(model.matrix(~ -1 + zhuClin$P_Stage))
rownames(zhuP) <- sub("zhuClin$", "", rownames(zhuP), fixed=T)
rownames(zhuP) <- sub(")", "_", rownames(zhuP), fixed=t)
z <- rbind(zhuExpr, zhuP)

# let pen be a vector of penalty factors for each coefficient. We want to preserve the pathological stage in the model, so we assign it 0 as penalty factor.
pen <- c(rep(1, nrow(dirExpr)), rep(0, nrow(dirP)))

# run the elastic net model
cv.fit <- cv.glmnet(x=t(x), y=factor(dirClin$os3yr), nfolds=10, alpha=.1, family="binomial", penalty.factor=pen)
plot(cv.fit)
fitEnet <- glmnet(x=t(x), y=factor(dirClin$os3yr), family="binomial", alpha=.1, lambda=cv.fit$lambda.min, penalty.factor=pen)

yhatEnet <- predict(fitEnet, t(z), type="response", s="lambda.min")
boxplot(yhatEnet ~ zhuClin$os3yr, ylab="3-year OS prediction (%)", xlab="3-year OS", main="elastic net - molecular + clinical features")
stripchart(yhatEnet ~ zhuClin$os3yr,pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

#  build a more robust classifier( bootstrapping)
tmp <- c()
try.cv.fit <- c()
tryfit <- c()
for(i in c(1:1000)) {
  print(i)
  N <- sample(colnames(dirExpr),274,replace=TRUE)
  try.cv.fit <- cv.glmnet(x=t(x[,N]), y=factor(dirClin[N,"os3yr"]), nfolds=10, alpha=.1, family="binomial", penalty.factor=pen)
  tryfit <- as.numeric(glmnet(x=t(x[,N]), y=factor(dirClin[N,"os3yr"]), family="binomial", alpha=.1, lambda=try.cv.fit$lambda.min, penalty.factor=pen)$beta)
  tmp <- cbind(tmp,tryfit)
  }

tmp <- apply(tmp,1,function(x){length(which(abs(x)>0))})
select_features <- which(tmp>=quantile(tmp,probs=.99))
rm(tmp1)
w <- as.data.frame(t(x[select_features,]))
w$os3yr <- dirClin$os3yr

# run the logit model on the top selected features
boostEnetfit <- glm(w$os3yr~. ,data = w, family = "binomial")
rm(w)

yhatboostEnet <- predict(object=boostEnetfit, newdata=as.data.frame(t(z)), type="response")
boxplot(yhatboostEnet ~ zhuClin$os3yr, ylab="prediction of 3-year OS probability (%)", xlab="3-year OS", main="elastic net - molecular + clinical features")
stripchart(yhatboostEnet ~ zhuClin$os3yr,pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

######################################################################################################################################
# 4. build the model based on clin + molecular features using RandomForest
######################################################################################################################################

# let's find the optimal value of mtry
a <- tuneRF(x=t(x), y=factor(dirClin$os3yr), stepFactor=1.5, doBest=TRUE, trace=TRUE, ntreeTry=1000,type="regression")

# run the model with this mtry value
fitRF <- randomForest(x=t(x), y=factor(dirClin$os3yr), mtry=a$mtry, do.trace=10, ntree=1000, importance=TRUE,type="regression")
yhatRF <- predict(fitRF, t(z), type="prob")[,2]

boxplot(yhatRF ~ zhuClin$os3yr, ylab="3-year OS prediction (%)", xlab="3-year OS", main= "Random Forest - molecular + clinical features")
stripchart(yhatRF ~ zhuClin$os3yr, pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

###################################################################################################################
# 5. principal component regression
##############################################################################################################

# make the formula Z~t(x) +pstage

fitPcr <- pcr(dirClin$os3yr ~ t(x), ncomp=10,validation = "CV", family="binomial")
yhatPcr <- predict(fitPcr, comps = 1:9,t(z), type="response")

boxplot(yhatPcr ~ zhuClin$os3yr, ylab="3-year OS prediction (%)", xlab="3-year OS", main= "Principal component regression - molecular + clinical features")
stripchart(yhatPcr ~ zhuClin$os3yr, pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

###################################################################################################################
# 6. partial leased square
##############################################################################################################
fitPls <- plsr(dirClin$os3yr ~ t(x), ncomp=10,validation = "CV", family="binomial")
yhatPls <- predict(fitPls, comps = 1:9,t(z), type="response")

boxplot(yhatPls ~ zhuClin$os3yr, ylab="3-year OS prediction (%)", xlab="3-year OS", main= "Principal component regression - molecular + clinical features")
stripchart(yhatPls ~ zhuClin$os3yr, pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

######################################################################################################################################
# 7. plot ROC curves to asses the performance of our models using pROC package
######################################################################################################################################

rocClin <- roc(predictor=as.numeric(yhatClin),response=as.numeric(zhuClin$os3yr),ci=TRUE)
rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.numeric(zhuClin$os3yr),ci=TRUE)
rocBoostEnet <- roc(predictor=as.numeric(yhatboostEnet), response=as.numeric(zhuClin$os3yr),ci=TRUE)
rocRF <- roc(predictor=as.numeric(yhatRF), response=as.numeric(zhuClin$os3yr),ci=TRUE)
rocPcr <- roc(predictor=as.numeric(yhatPcr), response=as.numeric(zhuClin$os3yr),ci=TRUE)
rocPls <- roc(predictor=as.numeric(yhatPls), response=as.numeric(zhuClin$os3yr),ci=TRUE)


plot.roc(rocClin,col="royalblue")
plot.roc(rocEnet,add=TRUE,col="red")
plot.roc(rocRF,add=TRUE,col="orange")
plot.roc(rocPcr,add=TRUE,col="darkgreen")
plot.roc(rocPls,add=TRUE,col="black")
plot.roc(rocBoostEnet,add=TRUE,col="green")

#concatenate AUC + 95CI 
txtClin <- paste("Logit Clin. AUC =",format(x=rocClin$auc, digits=2),", 95% CI:",format(x=as.numeric(rocClin$ci)[1],digits=2),"-",format(x=as.numeric(rocClin$ci)[3],digits=2))
txtEnet <- paste("Elastic Net AUC =",format(x=rocEnet$auc, digits=2),", 95% CI:",format(x=as.numeric(rocEnet$ci)[1],digits=2),"-",format(x=as.numeric(rocEnet$ci)[3],digits=2))
txtRF <- paste("Random Forest AUC =",format(x=rocRF$auc, digits=2),", 95% CI:",format(x=as.numeric(rocRF$ci)[1],digits=2),"-",format(x=as.numeric(rocRF$ci)[3],digits=2))
txtPcr <- paste("PrinComp AUC =",format(x=rocPcr$auc, digits=2),", 95% CI:",format(x=as.numeric(rocPcr$ci)[1],digits=2),"-",format(x=as.numeric(rocPcr$ci)[3],digits=2))
txtPls <- paste("Part. Least Sq. AUC =",format(x=rocPls$auc, digits=2),", 95% CI:",format(x=as.numeric(rocPls$ci)[1],digits=2),"-",format(x=as.numeric(rocPls$ci)[3],digits=2))
txtBoostEnet <- paste("Bootstrap Elastic Net AUC =",format(x=rocBoostEnet$auc, digits=2),", 95% CI:",format(x=as.numeric(rocBoostEnet$ci)[1],digits=2),"-",format(x=as.numeric(rocBoostEnet$ci)[3],digits=2))

title(main="performance of the models 
predicting the probability of 3 years OS",outer=TRUE)
text(x=.65, y=.25, labels=paste(txtClin), col="royalblue", adj=0)
text(x=.65, y=.2, labels=paste(txtEnet), col="red", adj=0)
text(x=.65, y=.15, labels=paste(txtRF), col="orange", adj=0)
text(x=.65, y=.1, labels=paste(txtPcr), col="darkgreen", adj=0)
text(x=.65, y=.05, labels=paste(txtPls), col="black", adj=0)
text(x=.65, y=.3, labels=paste(txtBoostEnet), col="green", adj=0)

# plot the roc curves separately
plot.roc(rocClin,col="royalblue", main="logit (clinical features only)")
plot.roc(rocEnet,col="red", main="Elastic Net (clinical + molecular features)")
plot.roc(rocRF,col="orange", main="Random Forest (clinical + molecular features)")
plot.roc(rocPcr,col="darkgreen", main="Principal Component Regression (clinical + molecular features)")
plot.roc(rocPls,col="black", main="Partial Least Square (clinical + molecular features)")
plot.roc(rocBoostEnet,col="green", main="Bootstrap Elastic Net (clinical + molecular features)")


######################################################################################################################################
# 8. draw the kaplan meier curves based on the predictors (high and low risk groups based on the median)
######################################################################################################################################

# for each model, let's assign to "high-risk" the group of patients whose yhat > median(yhat) 
# and to "low-risk" the group of patients whose yhat<median(yhat). We set high risk =1, low risk =0.

riskClin <- ifelse(yhatClin >= median(yhatClin), 1, 0)
riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))
names(riskEnet) <- rownames(yhatEnet)
riskboostEnet <- as.vector(ifelse(yhatboostEnet >= median(yhatboostEnet), 1, 0))
names(riskboostEnet) <- names(yhatboostEnet)
riskRF <- ifelse(yhatRF >= median(yhatRF), 1, 0)
riskPcr <- as.vector(ifelse(yhatPcr >= median(yhatPcr), 1, 0))
names(riskPcr) <- rownames(yhatPcr)
riskPls <- as.vector(ifelse(yhatPls >= median(yhatPls), 1, 0))
names(riskPls) <- rownames(yhatPls)

# for each model, draw the kaplan meir curves and compute the log rank test
par(mfrow=c(1,1))

plot(survfit(Surv(zhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhuClin$VITAL_STATUS) ~ riskClin), main="logit model", xlab="months",ylab="probability of OS",col= c("blue","magenta"), lwd=3)
survdiff(Surv(zhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhuClin$VITAL_STATUS) ~ riskClin, rho=0)
abline(v=36,col="red",lty=2,lwd=3)

plot(survfit(Surv(zhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhuClin$VITAL_STATUS) ~ riskEnet), main="elastic net model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(Surv(zhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhuClin$VITAL_STATUS) ~ riskEnet, rho=0)
abline(v=36,col="red",lty=2,lwd=3)

plot(survfit(Surv(zhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhuClin$VITAL_STATUS) ~ riskboostEnet), main=" Bootstrap elastic net model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(Surv(zhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhuClin$VITAL_STATUS) ~ riskboostEnet, rho=0)
abline(v=36,col="red",lty=2,lwd=3)

plot(survfit(Surv(zhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhuClin$VITAL_STATUS) ~ riskRF), main="random forest model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(Surv(zhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhuClin$VITAL_STATUS) ~ riskRF, rho=0)
abline(v=36,col="red",lty=2,lwd=3)

plot(survfit(Surv(zhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhuClin$VITAL_STATUS) ~ riskPcr), main="Prin. Comp. Regression model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(Surv(zhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhuClin$VITAL_STATUS) ~ riskPcr, rho=0)
abline(v=36,col="red",lty=2,lwd=3)

plot(survfit(Surv(zhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhuClin$VITAL_STATUS) ~ riskPls), main="Partial least square model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(Surv(zhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhuClin$VITAL_STATUS) ~ riskPls, rho=0)
abline(v=36,col="red",lty=2,lwd=3)

######################################################################################################################################
# 9. draw the heatmaps based on the predictors (high and low risk groups based on the median)
######################################################################################################################################
par(oma=c(4,2,4,4))

riskEnet1 <- sort(riskEnet)
heatmap.2(x=z[which(abs(fitEnet$beta)>0),names(riskEnet1)],
          trace="none",ColSideColors=c("blue","magenta")[riskEnet1+1],
          col=greenred(50),scale="row", 
          breaks = seq(-2.5, 2.5, len = 51),
          Colv=TRUE,main="Elastic Net")

riskboostEnet1 <- sort(riskboostEnet)
heatmap.2(x=z[which(abs(boostEnetfit$coefficients)>0),names(riskboostEnet1)],
          trace="none",ColSideColors=c("blue","magenta")[riskboostEnet1+1],
          col=greenred(50),scale="row", 
          breaks = seq(-2.5, 2.5, len = 51),
          Colv=TRUE,main="Bootstrap Elastic Net")

riskRF1 <- sort(riskRF)
heatmap.2(x=z[which(abs(importance(fitRF,type=1))>.18),
              names(riskRF1)],trace="none",ColSideColors=c("blue","magenta")[riskRF1+1],
          col=greenred(50),scale="row",
          breaks=seq(-2.5,2.5,len=51),
          Colv=TRUE, main="Random Forest")

riskPcr1 <- sort(riskPcr)
features_Pcr <- sort(abs(fitPcr$coefficients[,,1]),index.return=TRUE,decreasing=TRUE)$ix[1:100]
heatmap.2(x=z[features_Pcr,names(riskPcr1)],trace="none",
          ColSideColors=c("blue","magenta")[riskPcr1+1],col=greenred(50),
          scale="row",breaks=seq(-2.5,2.5,len=51),
          Colv=TRUE, main="Princ. Comp. Regression")

riskPls1 <- sort(riskPls)
features_Pls <- sort(abs(fitPls$coefficients[,,1]),index.return=TRUE,decreasing=TRUE)$ix[1:100]
heatmap.2(x=z[features_Pls,names(riskPls1)],trace="none",
          ColSideColors=c("blue","magenta")[riskPls1+1],
          col=greenred(50),scale="row",
          breaks=seq(-2.5,2.5,len=51),Colv=TRUE, main="Partial Least Square")



