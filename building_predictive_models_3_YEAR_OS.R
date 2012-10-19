# Charles Ferté,MD,MSc
# 10/11/2012
# Sage Bionetworks

# code to develop the model in the Director's dataset (TS) and validate it in the Zhu dataset (VS)
# we want to predict the 3 year propbability of OS

# load the packages that are neeeded
require(glmnet)
require(randomForest)
require(synapseClient)
require(caret)
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

# Probability of 3 year overall survival is our dependant variable (vectors of response): y_dir and y_zhu
zhu_clin$y_zhu <- y_zhu <-  ifelse(zhu_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36,1,0)
dir_clin$y_dir <-y_dir <-  ifelse(dir_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36,1,0)

# make data coherent for dir since there are two cel files that are not clinically annotated
tmp <- intersect(rownames(dir_clin),sampleNames(dir))
dir_clin <- dir_clin[tmp,]
dir <- exprs(dir)[,tmp]
rm(tmp)

# transform zhu eset in matrix
zhu <- exprs(zhu)

# get rid of the junk features
zhu <- zhu[-c(grep("AFFX",rownames(zhu))),]
dir <- dir[rownames(zhu),]

#check if there is any latent structure in the data (are the datasets comparable ?)
s <- svd(cbind(dir,zhu))
par(mfrow=c(2,2))
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
# 1. build a model based on clin variables of interest only (pStage, gender, age,smoking) (logisitic regression)
######################################################################################################################################
fit <- glm(dir_clin$y_dir~dir_clin$P_Stage + dir_clin$GENDER + dir_clin$Age,data=as.data.frame(t(dir_clin)),family="binomial")
summary(fit)

#predict in zhu
yhat <- predict(fit,newdata=as.data.frame(t(zhu_clin)),type="response")
yhat

## Brian there is a problem just above -> yhat should have only 62 values !!! and not 299 !
par(mfrow=c(1,1))
boxplot(yhat~zhu_clin$y_zhu,ylab="3-year OS prediction (%)",xlab="3-year OS",main="logit model - clinical features only")
stripchart(yhat~zhu_clin$y_zhu,pch=20,col="royalblue",vertical=TRUE,add=TRUE,cex=.6)

######################################################################################################################################
# 2. build a model based on logistic regression of clin + gene expression features
######################################################################################################################################

# first create datasets combining clin + molecular features in the same matrix 
# x to be the training set, z to be the validation set

x <- rbind(dir,as.numeric(as.factor(dir_clin$P_Stage)),as.numeric(as.factor(dir_clin$GENDER)),as.numeric(dir_clin$Age))
rownames(x)[(length(rownames(x))-2):length(rownames(x))] <- c("P_Stage","GENDER","Age")

z <- rbind(zhu,as.numeric(as.factor(zhu_clin$P_Stage)),as.numeric(as.factor(zhu_clin$GENDER)),as.numeric(zhu_clin$Age))
rownames(z)[(length(rownames(z))-2):length(rownames(z))] <- c("P_Stage","GENDER","Age")

fit <- glm(dir_clin$y_dir~.,data=as.data.frame(t(x)),family="binomial")

#predict in zhu
yhat <- predict(fit,newdata=as.data.frame(t(z)),type="response")
yhat

## Brian there is a problem just above -> yhat should have only 62 values !!! and not 299 !
par(mfrow=c(1,1))
boxplot(yhat~zhu_clin$y_zhu,ylab="3-year OS prediction (%)",xlab="3-year OS",main="logit model - clinical + molecular features")
stripchart(yhat~zhu_clin$y_zhu,pch=20,col="royalblue",vertical=TRUE,add=TRUE,cex=.6)

#########################################################################################################################################
# 3. build the model based on molecular features using elasticnet in a more ridge tend (alpha=.1) but not penalizing the clinical features
#########################################################################################################################################

set.seed(1234567)
pen <- c(rep(1,times=length(rownames(dir))),rep(0,times=3))
cv.fit <- cv.glmnet(x=t(x), y=factor(y_dir), nfolds=10, alpha=.1, family="binomial",penalty.factor=pen)
plot(cv.fit)
fit <- glmnet(x=t(x),y=factor(y_dir),family="binomial",alpha=.1,lambda=cv.fit$lambda.min,penalty.factor=pen)
table(as.numeric(fit$beta)!=0)    

yhat <- predict(fit,t(z),type="response",s="lambda.min")


par(mfrow=c(1,1))
boxplot(yhat~zhu_clin$y_zhu,ylab="3-year OS prediction (%)",xlab="3-year OS", main="elastic net")
stripchart(yhat~zhu_clin$y_zhu,pch=20,col="royalblue",vertical=TRUE,add=TRUE,cex=.6)

######################################################################################################################################
# 4. build the model based on clin + molecular features using RandomForest
######################################################################################################################################


fit <- randomForest(x=t(x),y=factor(y_dir),ntree=100, do.trace=10)
yhat <- predict(fit,t(z),type="prob")
boxplot(yhat~y_zhu,ylab="3-year OS prediction (%)",xlab="3-year OS",main= "Random Forest")
stripchart(yhat~y_zhu,pch=20,col="royalblue",vertical=TRUE,add=TRUE,cex=.6)

######################################################################################################################################
# plot a ROC curve to asses the performance of our model
######################################################################################################################################
require(ROCR)
Pred <- prediction(as.numeric(yhat),as.numeric(y_zhu))
Perf <- performance(prediction.obj=Pred,"tpr","fpr")
AUC <- performance(prediction.obj=Pred,"auc")
plot(Perf, col="royalblue",main="predicting KRAS G12C in ccle 
(modele trained in tcga + chemores)")
text(x=.7,y=.4,labels=paste("AUC=",format(x=AUC@y.values,digits=2)),col="royalblue")




