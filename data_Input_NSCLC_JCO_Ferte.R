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
