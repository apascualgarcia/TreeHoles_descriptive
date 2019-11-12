#################################################
#      SEM models with Lavaan                   #
#################################################
#
# This script considers a table of functions measured in [1]
# and allows performing Structural Equation Modelling analysis
# to study the relationships between these functions and the
# influence of the existence of community classes.
# The models used in the analysis presented in [2] are found in
# the folder data/LavaanModels and explained in the README file.
# The outputs of this script will be stored at data/SEM. The
# script will return the SEM obtained and a .txt file with the summary.
###############################################
# USAGE: Download the necessary data from 
# and locate it in the the folder data/source
# Create your SEM model and store it in data/LavaanModels or
# use some of the available ones. Edit your choices for
# the model below. 
################################################
# LICENSE: 
# This script is protected by the license provided in the repository:
# https://github.com/apascualgarcia/TreeHoles
################################################
# AP-G - apascualgarcia.github.io
# August 2019 Theoretical Biology, ETH-Zürich
################################################
# REFERENCES:
# [1] 
# [2] Pascual-García, A., & Bell, T. (2019). Community-level signatures of ecological succession 
#     in natural bacterial communities. bioRxiv, 636233.
################################################
#rm(list=ls())

#-- Load packages needed
library(gplots) 
library(MASS)
library(lavaan) #version 0.63 # 
library(plyr) # To use join_all

### START EDITING
# Set the working directory (path for the repo)
dirRoot="~/TreeHoles_descriptive/"
setwd(dirRoot)
# --- Parameters for the model:
rndMod=0 # If you want to run a model with the classes randomized fix to one
group=1 # Are you considering the different classes? if not=0
# ..... Fix the path and file
pathModel="data/LavaanModels/"
selectModel="Manifest_Milestone3.lav"
### STOP EDITING

# Read Files --------------------------------------------------------------
# .... The functions
pathSource="data/source/"
fileFun="20151016_Functions_remainder.csv"
# .... The properties of the functional groups you will use for your model 
fileClass="samples_metadata_time0.tsv"
# .... The output path (The files saved with the model are created according with "selectModel" variable)
pathOut="data/SEM"
#-- Read the clusters to which each sample is associated
setwd(pathSource)
metadata=read.table(fileClass,sep="\t")
Class=metadata[,c(1,4)] # samples and SparCC partition          
# ..... create a randomization of class
idxx=permute::shuffle(dim(Class)[1])
ClassRnd=Class[idxx,2]
Class=cbind(Class,ClassRnd)
colnames(Class)=c("Community","Id","Rnd")

#-- Read in functional data, and perform some operations
dd.data=read.csv(fileFun)
setwd(dirRoot)

#---- Exclude negative controls
dd.data=dd.data[dd.data$Community!="blank",]
dd.data=na.omit(dd.data)

#---- Take just  functions at time 7 and exclude normalized respiration
idxx=colnames(dd.data) 
dd.data=subset(dd.data,select=c(Community,Replicate,grep("7", idxx)))
dd.data=subset(dd.data,select= -c(pgRPC.7))

# Average and measurement error across replicas --------------------
#---- Average function measurements across the replicates
dd=aggregate(dd.data[3:ncol(dd.data)],list(dd.data$Community),function(x){mean(x,na.rm=T)})
row.names(dd)=unique(dd.data$Community)
colnames(dd)[1]="Community"

#---- log normalize
dd[,2:ncol(dd)]=log(dd[,2:ncol(dd)]+1) # Pseudocount

#---- Split function measurements across the replicates and look for common samples
dd.split=split(dd.data,dd.data$Replicate) # We first create a list split by replica
dd.join=join_all(dd.split,by="Community") # The problem is that not all the replicas have all samples, so first join them with a common community col
dd.join.clean=dd.join[complete.cases(dd.join),] # so we get just complete rows
idx.consens=dd.join.clean$Community # and now the names

Nrep=4
dd.split.log=list()
for(i in 1:Nrep){
  matched=match(idx.consens,dd.split[[i]]$Community)
  dd.split[[i]]=dd.split[[i]][matched,]
  #rownames(dd.split[[i]])=unique(dd.split[[i]]$Community)
  dd.split[[i]]=subset(dd.split[[i]],select=-c(Community,Replicate))
  dd.split.log[[i]]=log(dd.split[[i]]+1)
}

# ---- Compute measurement error (average of the all-against-all correlations)
Ncols=dim(dd.split.log[[1]])[2]
Npair=Nrep*(Nrep-1)/2
meas.err=matrix()
for(k in 1:Ncols){
  corr=0
  for(i in 1:(Nrep-1)){
    for(j in (i+1):Nrep){
      corr=corr+cor(dd.split.log[[i]][,k],dd.split.log[[j]][,k],use="complete.obs")
    }
  }
  corr=corr/Npair
  meas.err[k]=corr
}

#---- Build a new dataset with class 
dd.all=join_all(list(dd,Class),by="Community")
dd.all=dd.all[complete.cases(dd.all),]
dd.all$Id=factor(dd.all$Id)

# --- Build dummy columns from partition identities and an additional one with random identities
Id.f = dd.all$Id
dummies = model.matrix(~0+Id.f)
colnames(dummies)[1]="Id.f1"
Id.f.rnd= as.data.frame(sample(dd.all$Id)) # Generate random identities
colnames(Id.f.rnd)[1]="Id.rnd"

dd.all=cbind(dd.all,Id.f.rnd,dummies)

# Compute SEM --------------------------------
setwd(pathModel)
myModel <- readLines(selectModel)
setwd(dirRoot)
# ..... fit model
if(group==0){
  fit <- sem(myModel, data=dd.all)
}else{
  if(rndMod==0){ # Run a regular model
    selectModel=paste(selectModel,"group",sep=".") # Rename the label for the output files
    fit <- sem(myModel, data=dd.all, group="Id")
  }else{ # Run a randomization
    selectModel=paste(selectModel,"group","Rnd2",sep=".") # Rename the label for the output files
    fit <- sem(myModel, data=dd.all, group="Rnd")
  }

}

summary(fit, standardized=T,rsq=T, fit.measures = TRUE)
mod.fit=modindices(fit,sort. = TRUE)
#print(mod.fit[mod.fit$mi>3,])

# ...... save model
setwd(pathOut) # fix the output path
outModel=paste(selectModel,"mod",sep=".") # We save the model
save(myModel,file=outModel)

outFit=paste(selectModel,"fit",sep=".") # and the fit
sink(outFit) # Saving the summary is tricky, first I divert the output to the file
options(max.print=1000) # I skip the limitation in lines shown
summary(fit,standardized=T,rsq=T, fit.measures = TRUE) # I print the file
print(mod.fit)#[mod.fit$mi>3,]) # Also add modifications proposed
sink() # And close the pipe
setwd(dirRoot)
