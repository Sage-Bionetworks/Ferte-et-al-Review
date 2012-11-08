Ferté et al. J Clin Oncol. submission
================

Code used in Ferté et al., submitted to Journal of Clinical Oncology.

Herein, is described how readers that are computationally naive can actually interact with the data and reproduce the results displayed in this Review.

Objective:
The aim of this ancillary analysis is to address the specific bioinformatic issues of the development of molecular signatures as outlined in the text of this Review article.  The basis of this analysis is to utilize gene expression data from patients with non small-cell lung cancer (NSCLC) to predict the probability of 3 year overall survival in patients who underwent resection.

Selection of the datasets
To facilitate this analysis, gene expression datasets for patients diagnosed with early stage NSCLC were selected using the following criteria: each dataset had to be publicly available as raw data (CEL files and clinical summary), clinical covariates had to include pathological stage, adjuvant or neoadjuvant antineoplastic treatment status and overall survival information.  We selected and downladed the data of  four gene expression datasets that were publicly available:  Directors challenge data (available on caArray), Zhu et al (GEO: GSE14814), Hou et al (GEO: GSE19188), and TCGA LUSC (TCGA). 


Clinical data curation:
To ensure uniformity across the different datasets, patients that received any adjuvant or neoadjuvant treatment, patients without any overall survival information, as well as patients with stage IV or IIIB were excluded from the analysis.

Data storage:
The raw data (CEL files clinically curated with their companion clinicopathological information) were uploaded to the Synapse public portal (https://synapse.sagebase.org) for convenience under the following entity: Ferte et al. Review (syn87682). Extensive information on Synapse can be found at: 
https://sagebionetworks.jira.com/wiki/display/SYND/

Data handling:
All these data were handled using the R statistical software (http://www.r-project.org/) and selected packages that can be downloaded from Bioconductor or from the Comprehensive R Archive Network (CRAN).  Further information regarding the use of the Synapse R Client is available on the synapse wiki : 
https://sagebionetworks.jira.com/wiki/display/SYNR/R+Synapse+Client+Vignette.

Pre-processing:
We preprocessed the data with the most frequently used unsupervised normalization methods: rma, gcrma, mas5.0, dCHIP, frma, snm.

Statistical analyses: 
For educational and simplicity purposes, we intended to use the Director’s challenge data (rma normalized) as the training set and the Zhu et al. dataset (rma normalized) as the validation set. Logistic regression, Elastic net and Random forest were performed to predict the overall survival variables.

Code repository:
All code to perform these analyses are available via github at: https://github.com/Sage-Bionetworks/JCO-Review-Paper.  To learn more about using github, please see https://help.github.com/. 
After forking the github project, you will find a few different directories:
dataProcessing, dataAnalysis, Figures.  In the head directory, you will also find this readme file as well as a file named gitignore.  There may be some other files placed in the directory, depending on what terminal software you are using to run R.  
Files used to demonstrate the differences between supervised and unsupervised analysis (as described in Part 1 of the manuscript) are available in “dataProcessing” directory.  There are two files there, Norm.R , scalingDirZhu.R. In dataAnalysis is located predictiveModels.R. 
Start with Norm.R.  Each of the steps are well annotated within the code. Before running any code, we recommend installing the dependent packages used in the analysis:
affy http://www.bioconductor.org/packages/release/bioc/html/affy.html
gcrma http://www.bioconductor.org/packages/release/bioc/html/gcrma.html
frma http://www.bioconductor.org/packages/release/bioc/html/frma.html
hgu133afrmavecs http://www.bioconductor.org/packages/release/data/annotation/html/hgu133afrmavecs.html
hgu133plus2frmavecs http://www.bioconductor.org/packages/release/data/annotation/html/hgu133plus2frmavecs.html
hthgu133afrmavecs: http://www.mnmccall.com/software - http://www.mnmccall.com/software
snm http://www.bioconductor.org/packages/2.12/bioc/html/snm.html
synapseClient 'http://depot.sagebase.org/CRAN.R'
dependent packages:
caret
http://caret.r-forge.r-project.org/Classification_and_Regression_Training.html

Stepping through this file, you will find calls to pull in affymetrix files through the bioconductor affy package for each of the four datasets, following rma normalization, gcrma normalization, dCHIP normalization, MAS5 normalization, frma and barcode normalizations.  The commented codeblock at the end of this file demonstrates how to store the normalized data back in Synapse (though this is not necessary to run for the remainder of the code). The code then performs snm supervised normalization.  Comments at the end of each file provide links to store the models to Synapse.

Next, we go through the scalingDirZhu.R file to make the training and validation set comparable.

Finally we go through the predictiveModels.R file. The training and validation set (Dir and Zhu, respectively) used are rma normalized. Elastic Net, Random Forest, Logistic model, Partial Least square, principal component regression are implemented to predict 3 year overall survival.
The prediction performance is displayed through ROC AUC curves. Kaplan Meier and heatmaps are also drawn.
