########################################################
#   Community similarity within distance partitionings #
########################################################
#
# This is the sister script of LocationRandHierarchy2method.R 
# (please see that script for further details)
# in which we perform the same analysis creating first synthetic
# data to have a baseline against which compare the observed data.
# The default synthetic data consists of nested clusters of increasing size
# with decreasing similarity within the clusters. The mean and stdv
# of the distances generated within and between clusters can
# be controlled. The default nested
# design would be the one observed when clusters are defined
# joining samples that are geographically located at closer
# distances following for instance a random assemblage with
# some degree of dispersal limitation. 
# This nested structure have several variations in the way in which
# the mean and stdv of the distances are built, controlled by the variable
# setStructure, which can take the following values:
# a) "linear": the default structure
# b) "intermittent": Every two levels it decreases the mean to the minimum one.
# c) "random": After generating a linear structure, it shuffles the levels
# d) "inverse": change the order from maximum to minimum distances.
# e) "meanVar": Increase the stdv with increasing mean
# f) "stdvInverse": Decrease the stdv with increasing mean
################################################
# USAGE: Edit the parameters to create the matrix and select the
# desired method for the analysis. The main parameters are 
#
################################################
# The script generates a number of files and plots. The main files are ($label includes
# the method and distance used):
# 1) "DistanceDecay_ObsVsRandAndStdErr_$label" Has three columns, the value obtained
#    for the observed dataset with the selected method, and the mean value of the
#    randomizations for the immediately lower level, and the standard error.
# 2) "DistanceDecay_HierarchRandMean_$label" is a matrix with as many rows/cols as
#    spatial distances containing in the cell i,j the mean value of the selected
#    method across permutations, using the partitioning i when samples are permuted
#    within partitioning j. Therefore, it is an upper triangular matrix, but the 
#    observed values were added in the first column. A more convenient output for
#    plots is adding the observed values in the diagonal and using a lower triangular
#    matrix. This will be found in the file: 
#    -- "DistanceDecay_HierarchRandMean_$label.2ggplot.dat" which is also plotted in:
#    -- "DistanceDecay_HierarchRandMean_$label.pdf"
# 3) "DistanceDecay_HierarchRandStdErr_$label" a similar matrix as in 2) now containing
#    the standard errors. Since the observed values have no stdErr there is no input in
#    the first column.
# 4) If plotCtrl=1 and depending on the method used, a number of analysis will be plotted
#   in figures and the output of the analysis reported in files, including: anosim permutation
#   test analysis, analysis of the homogeneity of the variance across groups or Mantel tests.
################################################
# REFERENCES:
# [1] Pascual-García, A., & Bell, T. (2019). Community-level signatures of ecological succession 
#     in natural bacterial communities. bioRxiv, 636233.
# [2] Clarke, K. R. (1993). Non-parametric multivariate analysis of changes in 
#     community structure. Australian Journal of Ecology 18, 117–143.
# [3] Warton, D.I., Wright, T.W., Wang, Y. 2012. Distance-based multivariate analyses 
#     confound location and dispersion effects. Methods in Ecology and Evolution, 3, 89–101
################################################
# LICENSE: 
# This script is protected by the license provided in the repository:
# https://github.com/apascualgarcia/TreeHoles
################################################
# AP-G - apascualgarcia.github.io
# 1st version, June 25th, 2017 Silwood Park, Imperial College London
# 2nd version, August 1st, 2019 Theoretical Biology, ETH-Zürich
################################################

rm(list=ls())

library(vegan)  # For Mantel test and anosim
library(gplots)
library(ggplot2)
library(matrixStats) # for row/column stdv
library(reshape2)
library(gtools) # odd function

#**** START EDITING ****
# Set the working directory (path for the repo)
dirRoot="~/TreeHoles_descriptive/"
setwd(dirRoot)

# Setting parameters and input/output files ----------------------
# --- Parameters to create the random matrix
levels=8 # Number of levels in the interval [-1,1]. One gaussian distribution will be generated around each level
stdv=1 # Standard deviation
plotDist=1 # Should I create a heatmap with the generated matrix (=1)?

# --- Parameters for the anosim/adonis methods
Nperm0=9 # Number of permutations for the significance of the observed data
Nrnd=10 # Number of shuffled realizations for each nested level
Nperm1=1 # Number of permutations for the significance of the (already permuted) data. It may be interesting
# to use a large number to see if a permuted dataset is still significant. Otherwise, since we aim to test the
# significance of the observed data it is not retrieved, and can be low for computational efficiency.

# --- Output labels and directories 
#plotCtrl=1 # If =1 it will plot the permuted distro for each anosim analysis
distLabel="DistRand" # DistRand, since we are usingsynthetic distances
setLabel=paste("std",stdv,sep="") # Do not edit this
setMethod="adonis2" # anosim, mrpp or adonis2
setStructure="stdvInverse" # linear, intermittent, random, inverse, meanVar, stdvInverse

#**** STOP EDITING **** 

# Preliminaries ----------------
# --- Build labels and set working directory
if(setStructure=="meanVar"){
  stdv="increasing"
}else if(setStructure=="stdvInverse"){
  stdv="decreasing"
}
# ... build labels
labelOut=paste(distLabel,setLabel,setMethod,setStructure,sep="_")
dirOut="data/synthetic"

# The following would be used for the hierarchical rand procedure
fileMeanOut=paste("DistanceDecay_HierarchRandMean_",labelOut,".dat",sep="")
fileStdOut=paste("DistanceDecay_HierarchRandStdErr_",labelOut,".dat",sep="")
fileTableOut=paste("DistanceDecay_ObsVsRandAndStdErr_",labelOut,".dat",sep="")

setwd(dirOut)

# Create a random matrix ------------------
# --- Global values of the matrix
Nclus=2^seq(from=levels,to=1) # Number of clusters per level
Nclus=c(Nclus[1],Nclus) # the first level should be repeated to generate diagonal and off diagonal blocks
N=Nclus[1]*4 # total number of elements, four in each cluster in the finer classification
Nelem=N/Nclus # number of elements per cluster in each classification

# --- Fix the breaks along the interval, corresponding to the means of the gaussian distro
mean_int=2/(levels+1) # 2 corresponds to the size of the interval -1:1, where we will build our metric
stdv_int=mean_int/2
means=matrix(NA,1,(levels+1)) # this vector will contain the means of each level
stdv_vec=matrix(NA,1,(levels+1))
means[1]=1-mean_int/2 # from most similar
if((setStructure=="meanVar")||(setStructure == "stdvInverse")){
  stdv_vec[1]=stdv_int
}else{
  stdv_vec[1]=stdv
}
# --- Create the basic structure (linear model)
for(i in 2:(levels+1)){ # the breaks of the interval correspond to the mean of the correlations
  means[i]=means[i-1]-mean_int # decrease the similarity
  if((setStructure=="meanVar")||(setStructure == "stdvInverse")){
    stdv_vec[i]=stdv_vec[i-1]+stdv_int # increase the stdv constantly
  }else{
    stdv_vec[i]=stdv # keep it fixed
  }
}
# --- Modify it if other model is considered
if(setStructure=="intermittent"){ # every two levels come back to a high similarity
  for(i in 2:(levels+1)){ # the breaks of the interval correspond to the mean of the correlations
    if((odd(i)==TRUE)&& i != (levels+1)){
      means[i]=means[1]
    }
  }
}else if(setStructure == "random"){ # Randomize linear
  #setRand=sample.int((levels+1),(levels+1), replace = TRUE) # result 9 1 1 9 2 6 6 5 5
  setRand=c(9,1,1,9,2,6,6,5,5) # fixed for reproducibility
  means=means[setRand]
}else if(setStructure == "inverse"){ # invert the order
  idx=seq(from=(levels+1),to=1)
  means=means[idx]
}else if(setStructure == "stdvInverse"){ # invert for the stdv
  idx=seq(from=(levels+1),to=1)
  stdv_vec=stdv_vec[idx]
}
means # double ckeck
stdv_vec

# --- Build the correlation matrix
corrMat=matrix(NA,N,N) # Initialize correlation matrix
partMat=matrix(NA,N,(levels+1)) # Create a matrix of factors identifying the identity of the clusters
corrMatMean=matrix(0,(levels+1))
for(i in 1:(levels+1)){
  size=N/Nclus[i]
  if(i == 1){ # the first level has as many block distances as clusters
    Nblocks=Nclus[i]
  }else{ # and half off diagonal. The remainder clusters only have off diagonal blocks
    Nblocks=Nclus[i]/2
  }
  k=0
  Ncells=size*size
  for(j in 1:Nblocks){
    k=k+1 # this index will create a factor variable to identify the clusters
    M=matrix(rnorm(Ncells,mean=means[i],sd=stdv_vec[i]),size,size)
    M[which(M< -1)]=-1 # This can generate heavily U-shaped distros for high stdv
    M[which(M>1)]=1
    if(i == 1){ # The first level fills diagonal blocks
      idx=j-1
      col1=size*idx+1
      row1=col1
      clusInit=col1
    }else{ # all remaining blocks are off diagonal
      idx=2*j-1
      col1=size*idx+1
      row1=size*(idx-1)+1 
      clusInit=col1-size
    }
    col2=col1+size-1
    row2=row1+size-1
    clusEnd=col2
    corrMat[row1:row2,col1:col2]=M
    corrMat[col1:col2,row1:row2]=M
    partMat[clusInit:clusEnd,i]=k
    corrMatMean[i]=corrMatMean[i]+sum(corrMat[clusInit:clusEnd,clusInit:clusEnd])
  }
  if(i==1){
    norm=Ncells*Nblocks
  }else{
    norm=4*Ncells*Nblocks
  }
  corrMatMean[i]=corrMatMean[i]/norm
}
max(corrMat)  # check that the max is not larger than 1
min(corrMat) # check that the min is not lower than -1
DistCor=sqrt(2*(1-corrMat)) # Transform into distance
DistMeans=sqrt(2*(1-corrMatMean))
as.dist(DistCor)-> DistCor.dist
attr(DistCor.dist, "method") <- "dist"
setwd(dirRoot)

# --- Plot the distance matrix 
if(plotDist==1){
  setwd("/scripts")
  source("heatmap.2.mod.R")
  setwd(dirRoot)
  setwd(dirOut)
  figureOut=paste("heatmap_",labelOut,".pdf",sep="")
  pdf(figureOut,width=20,height = 20)
  heatmap.2.mod(DistCor[1:1024,1:1024],
                dendrogram = "none",scale="none",trace = "none",Rowv=NULL,Colv = NULL,
                xlab="Samples",ylab="Samples",
                cexRow = 0.01,cexCol=0.01,cex.lab=5,
                key.xlab = "Distance",key.title = "",
                key.par = list(cex.main=0.01, cex.lab=2.5,cex.axis=2,
                               mar=c(5,4.5,4,2)+0.1,mgp=c(3,1,0)), 
                #usr=c(-2,20,0,1)),xaxp=c(0,15,3),yaxp=c(0,0.15,3))ylog=TRUE,mgp=c(3,0.75,0),
                keysize=0.66,
                mgp=c(3,0.5,0),margins=c(5,5))
  dev.off()
}

# Compute ANOSIM, MRPP or ADONIS2 --------------------
anosimSumm=matrix(NA,nrow=levels,ncol=1)
anosimSig=matrix(NA,nrow=levels,ncol=1) # Will not be printed, since there is plotCtrl below for that
for(i in 1:levels){ # for each level
  if(setMethod=="anosim"){
    AA=anosim(DistCor.dist, as.factor(partMat[,i]), permutations = Nperm0)
    anosimSumm[i]=AA$statistic
    anosimSig[i]=AA$signif
  }else if(setMethod=="mrpp"){
    AA=mrpp(DistCor.dist, as.factor(partMat[,i]), permutations = Nperm0)
    anosimSumm[i]=AA$delta
    anosimSig[i]=AA$Pvalue
  }else{
    tmp.frame=as.data.frame(partMat[,i])
    colnames(tmp.frame)="X"
    AA=adonis2(DistCor.dist~X , tmp.frame, permutations = Nperm0)
    anosimSumm[i]=AA$F[1] 
    anosimSig[i]=AA$`Pr(>F)`[1]
  }
}
# ... Plot for the observed anosim
if(setMethod == "anosim"){
  methodLab="ANOSIM"
  ylabel=paste(methodLab," - R",sep="")
}else if(setMethod=="mrpp"){
  methodLab="MRPP"
  ylabel=paste(methodLab," - delta",sep="")
}else{
  methodLab="ADONIS"
  ylabel=paste(methodLab," - F",sep="")
}
df=data.frame(c(1:levels),DistMeans[1:levels],anosimSumm)
colnames(df)=c("Levels","Distance","Summary")

figureOut=paste("DistanceDecay-Obs_",labelOut,".pdf",sep="")
pdf(figureOut)
print(ggplot(df, aes(Distance,Summary))+
  geom_point()+geom_line()+
  ggtitle(paste("Stdv = ",stdv))+
  ylab(ylabel) +xlab("Distance Mean")+
  scale_size_manual(values=12)+
  theme_bw() +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 18),
        legend.title = element_text(size=18),
        legend.text=element_text(size=16),
        legend.key.size = unit(.5,"cm")))
dev.off()

figureOut=paste("DistanceDecay-Obs_byLevels_",labelOut,".pdf",sep="")
pdf(figureOut)
print(ggplot(df, aes(Levels,Summary))+
        geom_point()+geom_line()+
        ggtitle(paste("Stdv = ",stdv))+
        ylab(ylabel) +xlab("Levels")+
        scale_size_manual(values=12)+
        theme_bw() +
        theme(axis.title = element_text(size = 25),
              axis.text = element_text(size = 18),
              legend.title = element_text(size=18),
              legend.text=element_text(size=16),
              legend.key.size = unit(.5,"cm")))
dev.off()


# Shuffled ANOSIM across levels --------

# We store the observed values and we will overwrite with the mean randomizations below
AllSummaryMean=matrix(NA,nrow=(length(anosimSumm)+1),ncol=(levels+1)) # length(anosimSumm)=length(discuts)-1
AllSummaryMean[1:length(anosimSumm),1]=anosimSumm # Incorporate the observed values, leaving the last row NA
AllSummaryStd=matrix(NA,nrow=(length(anosimSumm)+1),ncol=(levels+1))
AllSummaryStd[1:length(anosimSumm),1]=matrix(0,nrow=length(anosimSumm),ncol=1) # There is no std for the observed anosim, so just equal to zero

# --- Now perform the hierarchical randomization
set.seed(1082018) # today
for(i in 1:(levels+1)){ # for each level
  if(i == 1){next} # We skip the first level (its constrained randomization is the full randomization), keeping the observed anosim in this col
  CTRL=how(blocks=partMat[,i]) # The shuffle will be made within blocks (see args of permute::shuffle)
  anosimSumm=matrix(NA, nrow=(i-1),ncol=Nrnd)
  anosimSig=matrix(NA, nrow=(i-1),ncol=Nrnd)
  for(u in 1:Nrnd){  # for each realization
    partMat.rnd=partMat
    idxx=permute::shuffle(dim(partMat)[1],control=CTRL)# it returns a list of shuffled indexes, we just need to know how many indexes
    partMat.rnd[,1:(i-1)]=partMat.rnd[idxx,1:(i-1)] # copy the partition with indexes shuffled
    for(k in 1:(i-1)){ # for each level above i  (shorter distances) perform ANOSIM permutting within i
      if(setMethod=="anosim"){
        AA=anosim(DistCor.dist, as.factor(partMat.rnd[,k]), permutations = Nperm1) # We do not care about the significance here
        anosimSumm[k,u]=AA$statistic
      }else if(setMethod=="mrpp"){
        AA=mrpp(DistCor.dist, as.factor(partMat.rnd[,k]), permutations = Nperm1)
        anosimSumm[k,u]=AA$delta
      }else{
        tmp.frame=as.data.frame(partMat.rnd[,k])
        colnames(tmp.frame)="X"
        AA=adonis2(DistCor.dist~X , tmp.frame, permutations = Nperm1)
        anosimSumm[k,u]=AA$F[1]
      }
    }
  }
  VecStd=rowSds(anosimSumm)/sqrt(Nrnd)
  VecMean=rowMeans(anosimSumm)
  # In the column i we have in rows the ANOSIM of the clusters of level i, randomizing the labels according to
  # the blocks found in the levels below i (labelled by the columns). In the first column we keep the observed
  # values, and therefore in row i there will be a jump with NA values for the columns k>i, since there is
  # no randomization made. See below where we transform the structure into data.frame after transposition, becoming more clear
  AllSummaryMean[1:(i-1),i]=VecMean
  AllSummaryStd[1:(i-1),i]=VecStd
}

# Create an output with the observed values, the random mean and stderrors
ObsVsRnd=AllSummaryMean[1:levels,1]
submatrix.tmp=AllSummaryMean[1:levels,2:dim(AllSummaryMean)[2]]
Rnd.tmp=diag(submatrix.tmp)#,AllSummaryMean[levels-1,1])
ObsVsRnd=cbind(ObsVsRnd,Rnd.tmp)
submatrix.tmp=AllSummaryStd[1:levels,2:dim(AllSummaryStd)[2]]
Rnd.tmp=diag(submatrix.tmp)#,AllSummaryStd[levels-1,1])
ObsVsRnd=cbind(DistMeans[1:levels],ObsVsRnd,Rnd.tmp) # The last mean in the vector
colnames(ObsVsRnd)=c("Mean distances",paste(setMethod,"Obs",sep=" "),paste(setMethod,"Rnd Mean",sep=" "),"StdErr")

# Print output -----------------------------------------
# Name files in the header of the program
# In the column i we have in rows the ANOSIM of the levels j<i (larger dists) keeping the clusters of i, except for the first column where we
# get the observed (not randomized values)
write.table(AllSummaryMean,sep="\t",file=fileMeanOut,quote=FALSE,row.names=FALSE)
write.table(AllSummaryStd,sep="\t",file=fileStdOut,quote=FALSE,row.names=TRUE)
write.table(ObsVsRnd,sep="\t",file=fileTableOut,quote=FALSE,row.names=TRUE)

# Figures ---------------------------------------------
# ... first some transformations to work with a df for ggplot
#AllSummaryMean=AllSummaryMean[1:10,1:11] # This is a temporal solution to a coding error above
# If you already run the script once and you want just to plot, load the data that will be stored below here
fileMeanOut2=paste("DistanceDecay_HierarchRandMean_",labelOut,".2ggplot.dat",sep="") # will create a file below

# .... Read from file (commented unless you read from file the results)
# file4ggplot=fileMeanOut2 # otherwise, read it if you generated in another run (this line and the 2 following should otherwise be commented)
# AllSummary=read.table(file=file4ggplot,header = TRUE)
# distcuts=seq(from=0.5,to=5.5,by=0.5) # Logarithmic scale base 10, spanning 5 orders of magnitude
# colnames(AllSummary)=c("distance",distcuts[1:levels-1]) # and names

# ... Create file
Obs=AllSummaryMean[,1] # observations are in the first colum
AllSummaryTmp=AllSummaryMean # we want them in the diagonal
AllSummaryTmp[,1]=NA # remove them
diag(AllSummaryTmp)=Obs # and add them in the diag
AllSummaryTmp=t(AllSummaryTmp) # and transpose


shapes=seq(1:levels)

colourSeries=c("black","blue","red","green4","grey",
               "cyan","orange","brown","violet","gold")

if(setStructure=="linear"){ # Only for the linear it makes sense to plot the means
  AllSummary=as.data.frame(cbind(DistMeans,AllSummaryTmp[,1:levels])) # add distances
  colnames(AllSummary)=c("distance",signif(DistMeans[1:levels],digits=2))
  xlabel="Mean distance within clusters"
  xvar="distance"
  df <- melt(AllSummary ,  id.vars = 'distance', variable.name = 'series')
  figureOut=paste("DistanceDecay-Shuffled_",labelOut,".pdf",sep="")
}else{
  AllSummary=as.data.frame(cbind(c(1:levels),AllSummaryTmp[,1:levels])) # add distances
  colnames(AllSummary)=c("levels",signif(DistMeans[1:levels],digits=2))
  xlabel="Levels"
  xvar="levels"
  df <- melt(AllSummary ,  id.vars = 'levels', variable.name = 'series')
  figureOut=paste("DistanceDecay-Shuffled_byLevels_",labelOut,".pdf",sep="")
}
yvar="value"
#seq(from=1,to=(levels))) # and names
pdf(file=figureOut, width=8)
print(ggplot(df, aes_string(xvar,yvar))+
  ylab(ylabel) +xlab(xlabel)+ggtitle(paste("Stdv = ",stdv))+
  geom_point(aes(shape=series,colour=series,fill=series),size=3)+
  geom_line(aes(colour = series))+
  scale_size_manual(values=12)+
  scale_shape_manual(values=shapes)+
  scale_colour_manual(values = colourSeries)+
  theme_bw() +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 18),
        legend.title = element_text(size=18),
        legend.text=element_text(size=16),
        legend.key.size = unit(.5,"cm")))
#  scale_color_brewer(palette="Paired")
dev.off()

# ... Save also this structure

write.table(AllSummary,sep="\t",file=fileMeanOut2,quote=FALSE,row.names=FALSE)
# 



