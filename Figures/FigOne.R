## Andrew Trister & Charles Ferté
## Sage Bionetworks
## Nov 15th 2012

## Raw Data PC plots - modified from Erich Huang's code
## FIRST: GENERATING UN-NORMALIZED, UN-BACKGROUND CORRECTED DATA
## We'll use the 'rma()' function with the normalize and background arguments set to false.

## REQUIRED LIBRARIES
require(Biobase)
require(affy)
require(snm)
require(ggplot2)
require(synapseClient)

synapseLogin()

## DEFINE getCelNames() FUNCTION
getCelNames <- function(x){
  require(affy)
  fileNames <- list.celfiles(path = x$cacheDir, full.names = TRUE)
}

## DEFINE getRawDat() FUNCTION
# Merely using the rma() function to extract summarized but not normalized or
# background-adjusted expression data
getRawDat <- function(x){
  require(affy)
  affyBatchObj <- ReadAffy(filenames = x)
  rawDat <- rma(affyBatchObj, normalize = FALSE, background = FALSE)
}

# Since these platforms are on two different Affymetrix platforms, we use
# the intersection of the Affymetrix probesets for comparative purposes

## DEFINE intersectFeatures() FUNCTION
## Take the entity list, get the feature names for the disparate Affy platforms and return
## the intersection of the feature names across platforms
intersectFeatures <- function(datMatList){
  featureList <- lapply(datMatList, rownames)
  intersectNames <- Reduce(intersect, featureList)
  return(intersectNames)
}


scaleFeatures <- function(datMatList){
  featuresList <- sapply(datMatList, function(x){
    t(scale(t(x),center=TRUE,scale=TRUE))
  })
  return(data.frame(featuresList))
}

# 
# normalize_to_X <- function(mean.x, sd.x, Y){
#  m.y <- rowMeans(Y)
#  sd.y <- apply(Y, 1, sd)
#  Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x 
#  Y.adj[sd.y == 0] <- mean.x[sd.y==0]
#  Y.adj
# }

# 
# scaleFeatures <- function(datMatList){
#   featuresList <- sapply(datMatList, function(x){
#     normalize_to_X(rowMeans(Zhu_snm),apply(Zhu_frma,1,sd),x)
#   })
#   return(data.frame(featuresList))
# }



#zhuExpr2 <- normalize_to_X(rowMeans(dirExpr), apply(dirExpr, 1, sd), zhuExpr)



## PULL IN THE RAW DATA FROM SYNAPSE
zhuRawEnt <- loadEntity('syn1439020')
houRawEnt <- loadEntity('syn1422295')
dirRawEnt <- loadEntity('syn1422422')
luscRawEnt <- loadEntity('syn1426948')


rawEntities <- list(zhuRawEnt, houRawEnt, dirRawEnt, luscRawEnt)
names(rawEntities) <- c('zhu', 'hou', 'dir', 'lusc')


celNamesList <- lapply(rawEntities, getCelNames)
rawDatList <- lapply(celNamesList, getRawDat)
rawDatMatList <- lapply(rawDatList, exprs)

commonFeatures <- intersectFeatures(rawDatMatList)

fullRawMat <- cbind(rawDatMatList$zhu[commonFeatures, ],
                    rawDatMatList$hou[commonFeatures, ],
                    rawDatMatList$dir[commonFeatures, ],
                    rawDatMatList$lusc[commonFeatures, ])

studyIndicator <- c(rep('zhu', ncol(rawDatMatList$zhu)),
                    rep('hou', ncol(rawDatMatList$hou)),
                    rep('dir', ncol(rawDatMatList$dir)),
                    rep('lusc', ncol(rawDatMatList$lusc)))

#Plot the principal components

svdObj <- svd(fullRawMat)
plot(svdObj$v[,1],svdObj$v[,2],
     col=c("royalblue","orange","aquamarine4","brown2")[as.factor(studyIndicator)],
     bg=c("royalblue","orange","aquamarine4","brown2")[as.factor(studyIndicator)],
     pch=c(8,24,23,19)[as.factor(studyIndicator)], font=2,
     cex=0.9, xlab="Principal component 1", ylab="Principal component 2")





# Now load the RMA normalized
zhuRMAEnt <- loadEntity('syn1436971')
Zhu_rma <- exprs(zhuRMAEnt$objects$Zhu_rma)

houRMAEnt <- loadEntity('syn1437174')
Hou_rma <- exprs(houRMAEnt$objects$Hou_rma)

dirRMAEnt <- loadEntity('syn1440819')
Dir_rma <- exprs(dirRMAEnt$objects$Dir_rma)

luscRMAEnt <- loadEntity('syn1457858')
Lusc_rma <- exprs(luscRMAEnt$objects$Lusc_rma)

RMADatMatList <- list(zhu =Zhu_rma,hou=Hou_rma,dir=Dir_rma,lusc=Lusc_rma)

RMAcommonFeatures <- intersectFeatures(RMADatMatList)

fullRMAMat <- cbind(RMADatMatList$zhu[RMAcommonFeatures, ],
                    RMADatMatList$hou[RMAcommonFeatures, ],
                    RMADatMatList$dir[RMAcommonFeatures, ],
                    RMADatMatList$lusc[RMAcommonFeatures, ])

RMAstudyIndicator <- c(rep('zhu', ncol(RMADatMatList$zhu)),
                    rep('hou', ncol(RMADatMatList$hou)),
                    rep('dir', ncol(RMADatMatList$dir)),
                    rep('lusc', ncol(RMADatMatList$lusc)))

#Plot the principal components

RMAsvdObj <- svd(fullRMAMat)
plot(RMAsvdObj$v[,1],RMAsvdObj$v[,2],
     col=c("royalblue","orange","aquamarine4","brown2")[as.factor(studyIndicator)],
     bg=c("royalblue","orange","aquamarine4","brown2")[as.factor(studyIndicator)],
     pch=c(8,24,23,19)[as.factor(studyIndicator)], font=2,
     cex=0.9, xlab="Principal component 1", ylab="Principal component 2")






# Now load the fRMA normalized
zhufRMAEnt <- loadEntity('syn1437066')
Zhu_frma <- exprs(zhufRMAEnt$objects$Zhu_frma)

houfRMAEnt <- loadEntity('syn1437182')
Hou_frma <- exprs(houfRMAEnt$objects$Hou_frma)

dirfRMAEnt <- loadEntity('syn1437194')
Dir_frma <- exprs(dirfRMAEnt$objects$Dir_frma)

luscfRMAEnt <- loadEntity('syn1458114')
Lusc_frma <- exprs(luscfRMAEnt$objects$Lusc_frma)

fRMADatMatList <- list(zhu =Zhu_frma,hou=Hou_frma,dir=Dir_frma,lusc=Lusc_frma)

fRMAcommonFeatures <- intersectFeatures(fRMADatMatList)

fullfRMAMat <- cbind(fRMADatMatList$zhu[fRMAcommonFeatures, ],
                    fRMADatMatList$hou[fRMAcommonFeatures, ],
                    fRMADatMatList$dir[fRMAcommonFeatures, ],
                    fRMADatMatList$lusc[fRMAcommonFeatures, ])

fRMAstudyIndicator <- c(rep('zhu', ncol(fRMADatMatList$zhu)),
                       rep('hou', ncol(fRMADatMatList$hou)),
                       rep('dir', ncol(fRMADatMatList$dir)),
                       rep('lusc', ncol(fRMADatMatList$lusc)))

#Plot the principal components

fRMAsvdObj <- svd(fullfRMAMat)
plot(fRMAsvdObj$v[,1],fRMAsvdObj$v[,2],
     col=c("royalblue","orange","aquamarine4","brown2")[as.factor(studyIndicator)],
     bg=c("royalblue","orange","aquamarine4","brown2")[as.factor(studyIndicator)],
     pch=c(8,24,23,19)[as.factor(studyIndicator)], font=2,
     cex=0.9, xlab="Principal component 1", ylab="Principal component 2")










# Now load the gcRMA normalized
zhugcRMAEnt <- loadEntity('syn1437007')
Zhu_gcrma <- exprs(zhugcRMAEnt$objects$Zhu_gcrma)

hougcRMAEnt <- loadEntity('syn1437176')
Hou_gcrma <- exprs(hougcRMAEnt$objects$Hou_gcrma)

dirgcRMAEnt <- loadEntity('syn1437188')
Dir_gcrma <- exprs(dirgcRMAEnt$objects$Dir_gcrma)

luscgcRMAEnt <- loadEntity('syn1457867')
Lusc_gcrma <- exprs(luscgcRMAEnt$objects$Lusc_gcrma)

gcRMADatMatList <- list(zhu =Zhu_gcrma,hou=Hou_gcrma,dir=Dir_gcrma,lusc=Lusc_gcrma)

gcRMAcommonFeatures <- intersectFeatures(gcRMADatMatList)

fullgcRMAMat <- cbind(gcRMADatMatList$zhu[gcRMAcommonFeatures, ],
                     gcRMADatMatList$hou[gcRMAcommonFeatures, ],
                     gcRMADatMatList$dir[gcRMAcommonFeatures, ],
                     gcRMADatMatList$lusc[gcRMAcommonFeatures, ])

gcRMAstudyIndicator <- c(rep('zhu', ncol(gcRMADatMatList$zhu)),
                        rep('hou', ncol(gcRMADatMatList$hou)),
                        rep('dir', ncol(gcRMADatMatList$dir)),
                        rep('lusc', ncol(gcRMADatMatList$lusc)))

#Plot the principal components

gcRMAsvdObj <- svd(fullgcRMAMat)
plot(gcRMAsvdObj$v[,1],gcRMAsvdObj$v[,2],
     col=c("royalblue","orange","aquamarine4","brown2")[as.factor(studyIndicator)],
     bg=c("royalblue","orange","aquamarine4","brown2")[as.factor(studyIndicator)],
     pch=c(8,24,23,19)[as.factor(studyIndicator)], font=2,
     cex=0.9, xlab="Principal component 1", ylab="Principal component 2")










# Now load the MAS5.0 normalized
zhuMAS5Ent <- loadEntity('syn1437053')
Zhu_MAS5 <- exprs(zhuMAS5Ent$objects$Zhu_MAS5)

houMAS5Ent <- loadEntity('syn1437178')
Hou_MAS5 <- exprs(houMAS5Ent$objects$Hou_MAS5)

dirMAS5Ent <- loadEntity('syn1437190')
Dir_MAS5 <- exprs(dirMAS5Ent$objects$Dir_MAS5)

luscMAS5Ent <- loadEntity('syn1457875')
Lusc_MAS5 <- exprs(luscMAS5Ent$objects$Lusc_MAS5)

MAS5DatMatList <- list(zhu =Zhu_MAS5,hou=Hou_MAS5,dir=Dir_MAS5,lusc=Lusc_MAS5)

MAS5commonFeatures <- intersectFeatures(MAS5DatMatList)

fullMAS5Mat <- cbind(MAS5DatMatList$zhu[MAS5commonFeatures, ],
                      MAS5DatMatList$hou[MAS5commonFeatures, ],
                      MAS5DatMatList$dir[MAS5commonFeatures, ],
                      MAS5DatMatList$lusc[MAS5commonFeatures, ])

MAS5studyIndicator <- c(rep('zhu', ncol(MAS5DatMatList$zhu)),
                         rep('hou', ncol(MAS5DatMatList$hou)),
                         rep('dir', ncol(MAS5DatMatList$dir)),
                         rep('lusc', ncol(MAS5DatMatList$lusc)))

#Plot the principal components

MAS5svdObj <- svd(fullMAS5Mat)
plot(MAS5svdObj$v[,1],MAS5svdObj$v[,2],
     col=c("royalblue","orange","aquamarine4","brown2")[as.factor(studyIndicator)],
     bg=c("royalblue","orange","aquamarine4","brown2")[as.factor(studyIndicator)],
     pch=c(8,24,23,19)[as.factor(studyIndicator)], font=2,
     cex=0.9, xlab="Principal component 1", ylab="Principal component 2")





# Now load the dCHIP normalized
zhudCHIPEnt <- loadEntity('syn1437063')
Zhu_dCHIP <- exprs(zhudCHIPEnt$objects$Zhu_dCHIP)

houdCHIPEnt <- loadEntity('syn1437180')
Hou_dCHIP <- exprs(houdCHIPEnt$objects$Hou_dCHIP)

dirdCHIPEnt <- loadEntity('syn1437192')
Dir_dCHIP <- exprs(dirdCHIPEnt$objects$Dir_dCHIP)

luscdCHIPEnt <- loadEntity('syn1457884')
Lusc_dCHIP <- exprs(luscdCHIPEnt$objects$Lusc_dCHIP)

dCHIPDatMatList <- list(zhu =Zhu_dCHIP,hou=Hou_dCHIP,dir=Dir_dCHIP,lusc=Lusc_dCHIP)

dCHIPcommonFeatures <- intersectFeatures(dCHIPDatMatList)

fulldCHIPMat <- cbind(dCHIPDatMatList$zhu[dCHIPcommonFeatures, ],
                     dCHIPDatMatList$hou[dCHIPcommonFeatures, ],
                     dCHIPDatMatList$dir[dCHIPcommonFeatures, ],
                     dCHIPDatMatList$lusc[dCHIPcommonFeatures, ])

dCHIPstudyIndicator <- c(rep('zhu', ncol(dCHIPDatMatList$zhu)),
                        rep('hou', ncol(dCHIPDatMatList$hou)),
                        rep('dir', ncol(dCHIPDatMatList$dir)),
                        rep('lusc', ncol(dCHIPDatMatList$lusc)))

#Plot the principal components

dCHIPsvdObj <- svd(fulldCHIPMat)
plot(dCHIPsvdObj$v[,1],dCHIPsvdObj$v[,2],
     col=c("royalblue","orange","aquamarine4","brown2")[as.factor(studyIndicator)],
     bg=c("royalblue","orange","aquamarine4","brown2")[as.factor(studyIndicator)],
     pch=c(8,24,23,19)[as.factor(studyIndicator)], font=2,
     cex=0.9, xlab="Principal component 1", ylab="Principal component 2")







# Now load the SNM normalized
zhusnmEnt <- loadEntity('syn1457384')
Zhu_snm <- exprs(zhusnmEnt$objects$Zhu_snm)

housnmEnt <- loadEntity('syn1457386')
Hou_snm <- exprs(housnmEnt$objects$Hou_snm)

dirsnmEnt <- loadEntity('syn1457380')
Dir_snm <- exprs(dirsnmEnt$objects$Dir_snm)

luscsnmEnt <- loadEntity('syn1457382')
Lusc_snm <- exprs(luscsnmEnt$objects$Lusc_snm)

snmDatMatList <- list(zhu =Zhu_snm,hou=Hou_snm,dir=Dir_snm,lusc=Lusc_snm)

snmcommonFeatures <- intersectFeatures(snmDatMatList)

fullsnmMat <- cbind(snmDatMatList$zhu[snmcommonFeatures, ],
                      snmDatMatList$hou[snmcommonFeatures, ],
                      snmDatMatList$dir[snmcommonFeatures, ],
                      snmDatMatList$lusc[snmcommonFeatures, ])

snmstudyIndicator <- c(rep('zhu', ncol(snmDatMatList$zhu)),
                         rep('hou', ncol(snmDatMatList$hou)),
                         rep('dir', ncol(snmDatMatList$dir)),
                         rep('lusc', ncol(snmDatMatList$lusc)))

#Plot the principal components

snmsvdObj <- svd(fullsnmMat)
plot(snmsvdObj$v[,1],snmsvdObj$v[,2],
     col=c("royalblue","orange","aquamarine4","brown2")[as.factor(studyIndicator)],
     bg=c("royalblue","orange","aquamarine4","brown2")[as.factor(studyIndicator)],
     pch=c(8,24,23,19)[as.factor(studyIndicator)], font=2,
     cex=0.9, xlab="Principal component 1", ylab="Principal component 2")





#plot the scaled SNM SVD


scaledsnmMat <-  scaleFeatures(list(a=snmDatMatList$zhu[snmcommonFeatures, ],
                                    b=snmDatMatList$dir[snmcommonFeatures, ]))

scaledsnmsvdObj <- svd(scaledsnmMat)
plot(scaledsnmsvdObj$v[,1],scaledsnmsvdObj$v[,2],
     col=c("royalblue","orange","aquamarine4","brown2")[as.factor(studyIndicator)],
     bg=c("royalblue","orange","aquamarine4","brown2")[as.factor(studyIndicator)],
     pch=c(8,24,23,19)[as.factor(studyIndicator)], font=2,
     cex=0.9, xlab="Principal component 1", ylab="Principal component 2")



#plot a table
plot(c(1,1,1,1),c(1,2,3,4),
     col=c("royalblue","orange","aquamarine4","brown2"),
     bg=c("royalblue","orange","aquamarine4","brown2"),
     pch=c(8,24,23,19), font=2,
     cex=0.9)


     
