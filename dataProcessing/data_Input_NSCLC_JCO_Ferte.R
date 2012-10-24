# Charles Ferté
# 09/07/2012
#Sage Bionetworks

# require packages
require(synapseClient)
synapseLogin(username="charles.ferte@sagebase.org",password="charles")

####################################################################################
# create a Zhu data entity inside the JCO Ferté project (syn87682)
####################################################################################
CELpath <- list.files("~/FELLOW/cferte/NSCLC_MA/CEL_BM/Zhu_CEL/", full.names=T,pattern=".CEL")
zhu_raw_exp <- Data(list(name = "Zhu CEL files", parentId = 'syn87682'))
zhu_raw_exp <- createEntity(zhu_raw_exp)

# add files into the data entity
zhu_raw_exp <- addFile(entity=zhu_raw_exp,CELpath)

# push the raw data into this entity
zhu_raw_exp <- storeEntity(entity=zhu_raw_exp)

####################################################################################
# create a Hou data entity inside the JCO Ferté project (syn87682)
####################################################################################
CELpath <- list.files("~/FELLOW/cferte/NSCLC_MA/CEL_BM/Hou_CEL/", full.names=T,pattern=".CEL")
Hou_raw_exp <- Data(list(name = "Hou CEL files", parentId = 'syn87682'))
Hou_raw_exp <- createEntity(Hou_raw_exp)

# add files into the data entity
Hou_raw_exp <- addFile(entity=Hou_raw_exp,CELpath)

# push the raw data into this entity
Hou_raw_exp <- storeEntity(entity=Hou_raw_exp)

####################################################################################
# create a Director's data entity inside the JCO Ferté project (syn87682)
####################################################################################
CELpath <- list.files("~/FELLOW/cferte/NSCLC_MA/CEL_BM/Dir_CEL/", full.names=T,pattern=".CEL")
Dir_raw_exp <- Data(list(name = "Directors CEL files", parentId = 'syn87682'))
Dir_raw_exp <- createEntity(Dir_raw_exp)

# add files into the data entity
Dir_raw_exp <- addFile(entity=Dir_raw_exp,CELpath)

# push the raw data into this entity
Dir_raw_exp <- storeEntity(entity=Dir_raw_exp)

####################################################################################
# create a TCGA LUSC's data entity inside the JCO Ferté project (syn87682)
####################################################################################
CELpath <- c(list.files("~/FELLOW/cferte/NSCLC_MA/NSCLC_Prognostic_Signature/LUSC_TCGA/CEL_files/broad.mit.edu_LUSC.HT_HG-U133A.Level_1.23.1005.0/", full.names=T,pattern=".CEL"),
             list.files("~/FELLOW/cferte/NSCLC_MA/NSCLC_Prognostic_Signature/LUSC_TCGA/CEL_files/broad.mit.edu_LUSC.HT_HG-U133A.Level_1.31.1005.0/", full.names=T,pattern=".CEL"),
             list.files("~/FELLOW/cferte/NSCLC_MA/NSCLC_Prognostic_Signature/LUSC_TCGA/CEL_files/broad.mit.edu_LUSC.HT_HG-U133A.Level_1.39.1005.0/", full.names=T,pattern=".CEL"),
             list.files("~/FELLOW/cferte/NSCLC_MA/NSCLC_Prognostic_Signature/LUSC_TCGA/CEL_files/broad.mit.edu_LUSC.HT_HG-U133A.Level_1.53.1005.0/", full.names=T,pattern=".CEL"),
             list.files("~/FELLOW/cferte/NSCLC_MA/NSCLC_Prognostic_Signature/LUSC_TCGA/CEL_files/broad.mit.edu_LUSC.HT_HG-U133A.Level_1.60.1005.0/", full.names=T,pattern=".CEL"))
             
LUSC_TCGA_raw_exp <- Data(list(name = "TCGA LUSC CEL files", parentId = 'syn87682'))
LUSC_TCGA_raw_exp <- createEntity(LUSC_TCGA_raw_exp)

# add files into the data entity
LUSC_TCGA_raw_exp <- addFile(entity=LUSC_TCGA_raw_exp,CELpath)

# push the raw data into this entity
LUSC_TCGA_raw_exp <- storeEntity(entity=LUSC_TCGA_raw_exp)

####################################################################################
# upload the clinical data files inside the JCO Ferté project (syn87682)
####################################################################################

load("~/FELLOW/cferte/NSCLC_MA/CLIN_DATA_FILES/TS_CLIN.Rdata")
load("~/FELLOW/cferte/NSCLC_MA/CLIN_DATA_FILES/VS_CLIN.Rdata")
load("~/FELLOW/cferte/NSCLC_MA/CLIN_DATA_FILES/VS2_CLIN.Rdata")
load("~/FELLOW/cferte/NSCLC_MA/CLIN_DATA_FILES/VS4_CLIN.Rdata")

###############################################################################

Dir_clin <- Data(list(name = "Directors clin files", parentId = 'syn87682'))
Dir_clin <- createEntity(Dir_clin)

# add files into the data entity
Dir_clin <- addObject(Dir_clin, DirClinF)

# push the raw data into this entity
Dir_clin <- storeEntity(entity=Dir_clin)

###############################################################################

Zhu_clin <- Data(list(name = "Zhu clin files", parentId = 'syn87682'))
Zhu_clin <- createEntity(Zhu_clin)

# add files into the data entity
Zhu_clin <- addObject(Zhu_clin,ZhuClinF)

# push the raw data into this entity
Zhu_clin <- storeEntity(entity=Zhu_clin)

###############################################################################

Hou_clin <- Data(list(name = "Hou clin files", parentId = 'syn87682'))
Hou_clin <- createEntity(Hou_clin)

# add files into the data entity
Hou_clin <- addObject(Hou_clin,HouClinF)

# push the raw data into this entity
Hou_clin <- storeEntity(entity=Hou_clin)

###############################################################################

Lusc_clin <- Data(list(name = "Lusc clin files", parentId = 'syn87682'))
Lusc_clin <- createEntity(Lusc_clin)

# add files into the data entity
Lusc_clin <- addObject(Lusc_clin,LuscClinF)

# push the raw data into this entity
Lusc_clin <- storeEntity(entity=Lusc_clin)

