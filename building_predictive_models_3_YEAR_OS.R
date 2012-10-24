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
synapseLogin(username="charles.ferte@sagebase.org",password="charles")


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
zhu_clin$y_zhu <-  ifelse(zhu_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36,1,ifelse(zhu_clin$VITAL_STATUS==0,NA,0))
zhu_clin <- zhu_clin[!is.na(zhu_clin$y_zhu),]
y_zhu  <- zhu_clin$y_zhu

dir_clin$y_dir <-  ifelse(dir_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36,1,ifelse(dir_clin$VITAL_STATUS==0,NA,0))
dir_clin <- dir_clin[!is.na(dir_clin$y_dir),]
y_dir <- dir_clin$y_dir

# make data coherent for dir since there are two cel files that are not clinically annotated
tmp <- intersect(rownames(dir_clin),sampleNames(dir))
dir_clin <- dir_clin[tmp,]
dir <- exprs(dir)[,tmp]
rm(tmp)

tmp <- intersect(rownames(zhu_clin),sampleNames(zhu))
zhu_clin <- zhu_clin[tmp,]
zhu <- exprs(zhu)[,tmp]
rm(tmp)

# get rid of the junk features
zhu <- zhu[-c(grep("AFFX",rownames(zhu))),]
dir <- dir[rownames(zhu),]

#check if there is any latent structure in the data (are the datasets comparable ?)
s <- svd(cbind(dir,zhu))

plot(s$v[,1],s$v[,2],col=c(rep("royalblue",times=299),rep("orange",times=62)),pch=20,xlab="PC1",ylab="PC2",main="initial datasets",cex=.5)

# quantile normalize zhu to have the same disribution than dir (make them comparable)
zhu <- normalize2Reference(zhu,rowMeans(dir))

#check if there is any latent structure in the data (are the datasets comparable ?)
s <- svd(cbind(dir,zhu))
plot(s$v[,1],s$v[,2],col=c(rep("royalblue",times=299),rep("orange",times=62)),pch=20,xlab="PC1",ylab="PC2",main="after quantile normalization",cex=.5)

# focus on the most variant probes
prob.var1  <- apply(zhu,1,var)
prob.var2 <- apply(dir,1,var)
mean.prob.var <- apply(cbind(prob.var1,prob.var2),1,mean)
tmp <- which(mean.prob.var>quantile(mean.prob.var,probs=.9))
dir <- dir[tmp,]
zhu <- zhu[tmp,]
rm(tmp,mean.prob.var,prob.var1,prob.var2)

#check if there is any latent structure in the data (are the datasets comparable ?)
s <- svd(cbind(dir,zhu))
plot(s$v[,1],s$v[,2],col=c(rep("royalblue",times=299),rep("orange",times=62)),pch=20,xlab="PC1",ylab="PC2",main="focus on the 10% top most variant features",cex=.5)

# make the two datasets have the same mean and variance
# Justin's function to rescale the VS to get the same mean/var than the TS
normalize_to_X <- function(mean.x, sd.x, Y){
  m.y <- rowMeans(Y)
  sd.y <- apply(Y, 1, sd)
  Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x 
  Y.adj[sd.y == 0] <- mean.x[sd.y==0]
  Y.adj
}

zhu <- normalize_to_X(rowMeans(dir),apply(dir,1,sd),zhu)

#check if there is any latent structure in the data (are the datasets comparable ?)
s <- svd(cbind(dir,zhu))
plot(s$v[,1],s$v[,2],col=c(rep("royalblue",times=299),rep("orange",times=62)),pch=20,xlab="PC1",ylab="PC2",main="align mean and variance of the two datasets",cex=.5)
rm(s)


######################################################################################################################################
# 2. build a model based on clin variables of interest only (pStage, gender, age) (logisitic regression)
######################################################################################################################################

dir_clin$P_Stage <- factor(dir_clin$P_Stage)
zhu_clin$P_Stage <- factor(zhu_clin$P_Stage, levels=levels(dir_clin$P_Stage))
dim(dir_clin)

fit <- glm(y_dir ~ P_Stage , data = dir_clin, family = "binomial")
summary(fit)
yhat <- predict(fit, zhu_clin, type = "response")

boxplot(yhat~y_zhu,ylab="3-year OS prediction (%)",xlab="3-year OS",main="logit model - clinical features only")
stripchart(yhat~y_zhu,pch=20,col="royalblue",cex=.6,vertical=TRUE,add=TRUE)

yhat1 <- yhat

#########################################################################################################################################
# 3. build the model based on molecular features using elasticnet 
#  not penalizing the clinical features in a more ridge setting (alpha=.1) 
#########################################################################################################################################

# let x and z be the datsets combining molecular features and P Stage
x <- rbind(dir,as.numeric(as.factor(dir_clin$P_Stage)))
rownames(x)[length(rownames(x))] <- "P_Stage"

z <- rbind(zhu,as.numeric(as.factor(zhu_clin$P_Stage)))
rownames(z)[length(rownames(z))] <- "P_Stage"

# let pen be a vector of penalty factors for each coefficient. We want to preserve the pathological stage in the model, so we assign it 0 as penalty factor.
pen <- c(rep(1,times=length(rownames(dir))),0)

# run the elastic net model
set.seed(1234567)
cv.fit <- cv.glmnet(x=t(x), y=factor(y_dir), nfolds=10, alpha=.1, family="binomial",penalty.factor=pen)
plot(cv.fit)
fit <- glmnet(x=t(x),y=factor(y_dir),family="binomial",alpha=.1,lambda=cv.fit$lambda.min,penalty.factor=pen)
yhat <- predict(fit,t(z),type="response",s="lambda.min")

boxplot(yhat~zhu_clin$y_zhu,ylab="3-year OS prediction (%)",xlab="3-year OS", main="elastic net - molecular + clinical features")
stripchart(yhat~zhu_clin$y_zhu,pch=20,col="royalblue",vertical=TRUE,add=TRUE,cex=.6)

yhat2 <- yhat

######################################################################################################################################
# 4. build the model based on clin + molecular features using RandomForest
######################################################################################################################################

# let's find the optimal value of mtry
a <- tuneRF(x=t(x),y=factor(y_dir),stepFactor=1.5,doBest=TRUE,trace=TRUE,ntreeTry=1000)

# run the model with this mtry value
set.seed(1234567)
fit <- randomForest(x=t(x),y=factor(y_dir),mtry=a$mtry, do.trace=10,ntree=1000, importance=TRUE)
yhat <- predict(fit,t(z),type="prob")[,2]

boxplot(yhat~y_zhu,ylab="3-year OS prediction (%)",xlab="3-year OS",main= "Random Forest - molecular + clinical features")
stripchart(yhat~y_zhu,pch=20,col="royalblue",vertical=TRUE,add=TRUE,cex=.6)

yhat3 <- yhat

######################################################################################################################################
# 5. plot a ROC curve to asses the performance of our models
######################################################################################################################################

require(ROCR)

Pred <- prediction(as.numeric(yhat1),as.numeric(y_zhu))
Perf <- performance(prediction.obj=Pred,"tpr","fpr")
AUC <- performance(prediction.obj=Pred,"auc")
plot(Perf, col="royalblue",main="performance of the models",lwd=2)
text(x=.35,y=.3,labels=paste("AUC logit clin. =",format(x=AUC@y.values,digits=2)),col="royalblue",adj=0,cex=.8)

Pred <- prediction(as.numeric(yhat2),as.numeric(y_zhu))
Perf <- performance(prediction.obj=Pred,"tpr","fpr")
AUC <- performance(prediction.obj=Pred,"auc")
plot(Perf, col="orange",lwd=2,add=TRUE)
text(x=.35,y=.25,labels=paste("AUC Elastic Net=",format(x=AUC@y.values,digits=2)),col="orange",adj=0,cex=.8)

Pred <- prediction(as.numeric(yhat3),as.numeric(y_zhu))
Perf <- performance(prediction.obj=Pred,"tpr","fpr")
AUC <- performance(prediction.obj=Pred,"auc")
plot(Perf, col="red",add=TRUE,lwd=2)
text(x=.35,y=.2,labels=paste("AUC Random Forest=",format(x=AUC@y.values,digits=2)),col="red",adj=0,cex=.8)

######################################################################################################################################
# 6. draw the kaplan meier curves based on the predictors (high and low risk groups based on the median)
######################################################################################################################################

# for each model, let's assign to "high-risk" the group of patients whose yhat > median(yhat) 
# and to "low-risk" the group of patients whose yhat<median(yhat). We set high risk =1, low risk =0.
risk1 <- ifelse(yhat1>=median(yhat1),1,0)
risk2 <- ifelse(yhat2>=median(yhat2),1,0)
risk3 <- ifelse(yhat3>=median(yhat3),1,0)

plot(survfit(Surv(zhu_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhu_clin$VITAL_STATUS)~risk1,data=zhu_clin), main="logit model: clinical variables only")
survdiff(Surv(zhu_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhu_clin$VITAL_STATUS)~risk1,data=zhu_clin,rho=0)

plot(survfit(Surv(zhu_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhu_clin$VITAL_STATUS)~risk2,data=zhu_clin),main="elastic net model: clinical + molecular features")
survdiff(Surv(zhu_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhu_clin$VITAL_STATUS)~risk2,data=zhu_clin,rho=0)

plot(survfit(Surv(zhu_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhu_clin$VITAL_STATUS)~risk3,data=zhu_clin), main="random forest model: clinical + molecular features")
survdiff(Surv(zhu_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH,zhu_clin$VITAL_STATUS)~risk3,data=zhu_clin,rho=0)