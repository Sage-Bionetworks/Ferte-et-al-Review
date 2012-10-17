# Charles Ferté
# 10/11/2012
#Sage Bionetworks

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

# make data coherent for dir since there are two cel files that are not clinically annotated
tmp <- intersect(rownames(dir_clin),sampleNames(dir))
dir_clin <- dir_clin[tmp,]
dir <- dir[,tmp]
rm(tmp)

# quantile normalize zhu to have the same disribution than dir (make them comparable)
zhu <- normalize2Reference(exprs(zhu),rowMeans(exprs(dir)))

# create vectors of response( y_dir and y_zhu) of 3 year overall survival (1 alive, 0 deceased)
zhu_clin$y_zhu <- y_zhu <-  ifelse(zhu_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36,1,0)
dir_clin$y_dir <-y_dir <-  ifelse(dir_clin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36,1,0)

######################################################################################################################################
# first build the model based on clin variables of interest only (pStage, gender, age,smoking) (linear regression)
######################################################################################################################################

fit <- glm(dir_clin$y_dir~dir_clin$P_Stage + dir_clin$GENDER + dir_clin$Age,data=as.data.frame(t(dir_clin)),family="binomial")
summary(fit)
confint.default(fit)
boxplot(fit$fitted.values~dir_clin$y_dir,ylab="3-year OS prediction (%)",xlab="3-year OS")
stripchart(fit$fitted.values~dir_clin$y_dir,pch=20,col="royalblue",vertical=TRUE,add=TRUE,cex=.6)

# predict in zhu
yhat_zhu <- predict.glm(fit,newdata=as.data.frame(zhu_clin),type="response",na.action = na.omit)

######################################################################################################################################
# second build the model based on linear regression of clin + ge features
######################################################################################################################################



######################################################################################################################################
# third build the model based on molecular features using elasticnet in a more ridge tend (alpha=.1)
######################################################################################################################################

# make zhu have the same mean and var that dir
set.seed(1234567)
identical(colnames(exprs(dir)),rownames(dir_clin))
x <- rbind(exprs(dir),as.numeric(as.factor(dir_clin$P_Stage)))
rownames(x)[22283] <- "P_Stage"
pen <- c(rep(1,times=length(rownames(exprs(dir)))),rep(0,times=1))

cv.fit <- cv.glmnet(x=t(exprs(dir)), y=factor(y_dir), nfolds=10, alpha=.1, family="binomial",penalty.factor=pen)
plot(cv.fit)
fit <- glmnet(x=t(exprs(dir)),y=factor(y_dir),family="binomial",alpha=.1,lambda=cv.fit$lambda.min,penalty.factor=pen)
table(as.numeric(fit$beta)!=0)    
yhat <- predict(fit, t(zhu),type="response")
boxplot(yhat~zhu_clin$y_zhu,ylab="3-year OS prediction (%)",xlab="3-year OS")
stripchart(yhat~zhu_clin$y_zhu,pch=20,col="royalblue",vertical=TRUE,add=TRUE,cex=.6)



# third build the model based on clin + molecular features using RandomForest

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


# adress the same results of AUC in boxplots with their IC95

