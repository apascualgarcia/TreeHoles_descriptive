########################################################
#   Community similarity within distance partitionings #
########################################################
#
# This script takes a distance matrix comparing the beta
# diversity of different samples, and another distance 
# matrix monitoring their geographic (Haversine) distance,
# and estimates the autocorrelation of the samples comparing
# the similarity of groups of samples closer than certain 
# spatial threshold, with respect to the between-groups 
# similarity. To estimate the significance of this similarity
# three methods are implemented, anosim, mrpp and adonis. 
# We aim to understand if the communities similarity changes across
# different distance thresholds. However, we deal with a nested design
# in which clusters at lower distances are nested within
# those at larger distances and, therefore, to understand
# which signal we should expect in a scenario of distance-decay
# (which would result for instance if neutral processes and
# dispersal limitation shape the communities) we also provide
# a sister script in which we use synthetic data with a decay
# in the dissimilarity of the communities at larger clustering
# thresholds. From that analysis we observed that anosim and mrpp
# capture the decay, anosim being more sensitive to an overlapping
# variance in the distances of the clusters (see [2]), while adonis
# just detect that the treatment (a.k.a the existence of clusters)
# plays a significant role, but it cannot detect the decay.
# In addition, in real data the size of the clusters is heterogeneous
# and also the variance of the distances within clusters, and
# to understand if this fact has any effect we studied how
# the observed value of the method used changes when the samples are
# permuted within clusters at larger distance thresholds
# (i.e. a hierarchical randomization). This allow us to control
# if any variation in the signal is due to the finer partitioning
# occuring at closer distances (more clusters), which will affect the
# assumption of the homogeneity of variance made by these
# methods, or due to a higher similarity of the samples (see [2], [3]). The rational 
# is that, if given a distance threshold, the similarity of the samples
# do not increase at lower distances, if we group the samples at
# lower distances the signal should not increase, and if an increase
# is observed it should be attributed to the new partitining. This
# kind of nested permutations can be done in adonis easily with 
# two levels, but if the number of levels increases the design 
# becomes very complicated so we perform it here manually. 
################################################
# USAGE: This script was designed for the dataset presented in [1], so
# you should edit the script to consider other matrices and thresholds 
# to be loaded at the distcut structure. The method is, otherwise, 
# compatible with other datsets. For the dataset in [1], you should
# edit just the type of matrix and method you want to use. You may
# want to edit in the code the number of randomizations or permutations. 
# The option plotCtrl will additionally generate a number of analysis
# such as plots of the permutation tests or analysis of the dispersion.
# For large numbers the method may take several hours in a regular
# computer.
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
library(stringi)
library(dendextend) # for dendrograms
library(permute)
library(ggplot2)
library(ade4)
library(matrixStats) # for row/column stdv
library(lattice) # to combine plots in the same...
library(gridExtra) # ... window with grid.arrange
library(reshape2)

## START EDITING
# Fix parameters ------------------
# Set the working directory (path for the repo)
dirRoot="~/TreeHoles_descriptive/"
setwd(dirRoot)
# Randomization parameters
Nperm0=999 # Number of permutations for the significance of the observed data
Nrnd=200 # Number of permutations for each nested level
Nperm1=1 # Number of permutations for the significance of the (already permuted) data. It may be interesting
         # to use a large number to see if a permuted dataset is still significant over the shufflings. Otherwise, since we aim to test the
         # significance of the observed data it is not retrieved, and can be low for computational efficiency.
# Fix some output labels and directories 
plotCtrl=1 # If =1 it will perform and plot additional tests, such as the permuted distro for each anosim analysis
distLabel="SparCC" # SparCC or SJD
setLabel="Time0" # Do not edit this
setMethod="anosim" # anosim, mrpp or adonis2
## STOP EDITING

# Load data --------------------------------------
# ... build labels
labelOut=paste(distLabel,setLabel,setMethod,sep="_")
fileMeanOut=paste("DistanceDecay_HierarchRandMean_",labelOut,".dat",sep="")
fileStdOut=paste("DistanceDecay_HierarchRandStdErr_",labelOut,".dat",sep="")
fileTableOut=paste("DistanceDecay_ObsVsRandAndStdErr_",labelOut,".dat",sep="")

# Load data and match samples  -----------------------------------------------
# --- Load GPS Haversine distances
dirSource="data/source"
setwd(dirSource)
inputGPS='Dist_GPS-Haversine.dat'
GPS.dist=read.table(inputGPS,sep=" ",header=TRUE, row.names=1)

# --- Load data for communities similarity and do some processing (SparCC matrix)
if(distLabel == "SparCC"){
   inputCor='corMat-SparCC_20151016_OTU_remainder.clean.samples.txt'
   CorSamples=read.table(inputCor,sep="\t",header=TRUE,row.names = 1)
}else{ # Jensen-Shannon
  inputCor='distMat_ShannonJensen_Samples_Time0.clean.dat'
  CorSamples=read.table(inputCor)
  DistCor=CorSamples # Already a distance
}

# --- Match distances
matched=match(rownames(CorSamples),rownames(GPS.dist)) # Look for the samples sequenced
matched=matched[!is.na(matched)] # delete those absent
GPS.dist=GPS.dist[matched,matched] # keep only those present in the first matrix
matched=match(rownames(CorSamples),rownames(GPS.dist)) # match again
CorSamples=CorSamples[!is.na(matched),!is.na(matched)] # and remove in the second matrix
matched=match(rownames(CorSamples),rownames(GPS.dist)) # final match, should not be needed but just in case
GPS.dist=GPS.dist[matched,matched] # we guarantee that the samples are in the same order in both matrices
# --- Once matched we transform into a distance
if(distLabel == "SparCC"){
  DistCor=sqrt(2*(1-CorSamples)) # If SparCC, convert into distances
}else{
  DistCor=CorSamples
}
as.dist(DistCor)-> DistCor.dist
attr(DistCor.dist, "method") <- "dist"
setwd(dirRoot)

# Clustering --------------------------------------------------------------
# We want to cluster samples and to cut the tree at different levels, to obtain different partitions by distance
distcuts=seq(from=0.5,to=5.5,by=0.5) # Logarithmic scale base 10, spanning 5 orders of magnitude 

# We now cluster samples 
clustGPS=hclust(as.dist(log10(GPS.dist+1)),method="single") # I use logarithmic distances to better represent it
dend1=as.dendrogram(clustGPS)

# and we cut the tree at different levels, to obtain different partitions
N=dim(GPS.dist)[1]
GPS.dist.cutree <- data.frame(rep(NA, N)) #, txt=rep("", N),  # as many cols as you need
maxPart=matrix()
i=0
for(hi in distcuts){
  i=i+1
  GPS.dist.cutree[,i]=cutree(dend1,h=hi)
  maxPart[i]=unique(GPS.dist.cutree[,i])[length(unique(GPS.dist.cutree[,i]))] # Max number of clusters
}
colnames(GPS.dist.cutree)=as.character(distcuts)
rownames(GPS.dist.cutree)=rownames(GPS.dist)

# ANOSIM -ADONIS analysis ------------------------------------------
# .... in the following the structures are called anosim even if adonis or mrpp could be used.
# --- First compute the observed anosim aross levels
dirOut="data/ANOSIM"
setwd(dirOut) # We start generating figures
anosimSumm=matrix(NA,nrow=(length(distcuts)-1),ncol=1)
anosimSig=matrix(NA,nrow=(length(distcuts)-1),ncol=1) # Will not be printed, since there is plotCtrl below for that
i=0
for(hi in distcuts){ # for each level
  i=i+1
  if(i==length(distcuts)){next}
  if(setMethod=="anosim"){
    AA=anosim(DistCor.dist, as.factor(GPS.dist.cutree[,i]), permutations = Nperm0)
    anosimSumm[i]=AA$statistic
    anosimSig[i]=AA$signif
  }else if(setMethod=="mrpp"){
    AA=mrpp(DistCor.dist, as.factor(GPS.dist.cutree[,i]), permutations = Nperm0)
    anosimSumm[i]=AA$delta
    anosimSig[i]=AA$Pvalue
  }else{
    tmp.frame=as.data.frame(GPS.dist.cutree[,i])
    colnames(tmp.frame)="X"
    AA=adonis2(DistCor.dist~X , tmp.frame, permutations = Nperm0)
    anosimSumm[i]=AA$F[1] 
    anosimSig[i]=AA$`Pr(>F)`[1]
  }
  if(i==length(distcuts)){next} # Only one cluster in the last step
  if(plotCtrl==1){ # Perform the standard controls for the observed values
    if(setMethod=="anosim"){ # controls specific of ANOSIM
      perm <- permustats(AA)
      outPlot=paste("anosimSignif_logDist-",hi,"_",distLabel,"_",setLabel,".pdf",sep="")
      pdf(outPlot,height=7,width=14)
      g1 = densityplot(perm)
      g2 = qqmath(perm)
      grid.arrange(g1,g2,nrow = 1,top=paste("distance = 10^",hi," m.",sep=""))
      dev.off()
      outFit=paste("anosimSummary_logDist-",hi,"_",distLabel,"_",setLabel,".dat",sep="") # and the fit
      sink(outFit) # divert the output to the file
      options(max.print=1000) #  skip the limitation in lines s
      print(summary(perm)) # print the file
      sink() # close the pipe
    }
    ## Calculate multivariate dispersions
    mod=betadisper(DistCor.dist,group = as.factor(GPS.dist.cutree[,i]))
    outFit=paste("betadisperSumm_logDist-",hi,"_",distLabel,"_",setLabel,".dat",sep="") # and the fit
    sink(outFit) # divert the output to the file
    print(anova(mod)) ## Perform test
    #print(permutest(mod, pairwise = TRUE, permutations = 999)) ## Permutation test for F
    print((mod.HSD <- TukeyHSD(mod)))     ## Tukey's Honest Significant Differences
    sink()
    outPlot=paste("betadisperSignif_logDis t-",hi,"_",distLabel,"_",setLabel,".pdf",sep="")
    pdf(outPlot,height=7,width=14)
    par(mfrow=c(1,2))
    plot(mod, ellipse = TRUE, hull = FALSE, main=distLabel,sub="")#,label.cex=0.1) # Plot the groups and distances to centroids on the first two PCoA axes, 1 sd data ellipse
    boxplot(mod,xlab="Community Class")
    dev.off()
  }
}
# We store the observed values and we will overwrite with the mean randomizations below
AllSummaryMean=matrix(NA,nrow=length(anosimSumm),ncol=length(distcuts)) # length(anosimSumm)=length(discuts)-1
AllSummaryMean[,1]=anosimSumm # Incorporate the observed values
AllSummaryStd=matrix(NA,nrow=length(anosimSumm),ncol=length(distcuts))
AllSummaryStd[,1]=matrix(0,nrow=length(anosimSumm),ncol=1) # There is no std for the observed anosim, so just equal to zero

# --- Prepare a figure for the observed values
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
df=data.frame(distcuts[1:(length(distcuts)-1)],anosimSumm)
colnames(df)=c("Distance","Summary")

figureOut=paste("DistanceDecay-Obs_",labelOut,".pdf",sep="")
pdf(figureOut)
print(ggplot(df, aes(Distance,Summary))+
  geom_point()+geom_line()+
  #ggtitle(paste("Stdv = ",stdv))+
  ylab(ylabel) +xlab("log10(distance) (meters)")+
  scale_size_manual(values=12)+
  theme_bw() +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 18),
        legend.title = element_text(size=18),
        legend.text=element_text(size=16),
        legend.key.size = unit(.5,"cm")))
dev.off()

# --- Now perform the hierarchical randomization
set.seed(1082018) # today
i=0
for(hi in distcuts){ # for each level
  i=i+1
  if(hi == 0.5){next} # We skip the first level (its constrained randomization is the full randomization), keeping the observed anosim in this col
  CTRL=how(blocks=GPS.dist.cutree[,i]) # The shuffle will be made within blocks (see args of permute::shuffle)
  anosimSumm=matrix(NA, nrow=(i-1),ncol=Nrnd)
  anosimSig=matrix(NA, nrow=(i-1),ncol=Nrnd)
  for(u in 1:Nrnd){  # for each realization
    GPS.dist.cutree.rnd=GPS.dist.cutree
    idxx=permute::shuffle(dim(GPS.dist.cutree)[1],control=CTRL)# it returns a list of shuffled indexes, we just need to know how many indexes
    GPS.dist.cutree.rnd[,1:(i-1)]=GPS.dist.cutree.rnd[idxx,1:(i-1)] # copy the partition with indexes shuffled
    for(k in 1:(i-1)){ # for each level above i  (shorter distances) perform ANOSIM permutting within i
      if(setMethod=="anosim"){
        AA=anosim(DistCor.dist, as.factor(GPS.dist.cutree.rnd[,k]), permutations = Nperm1) # We do not care about the significance here
        anosimSumm[k,u]=AA$statistic
      }else if(setMethod=="mrpp"){
        AA=mrpp(DistCor.dist, as.factor(GPS.dist.cutree.rnd[,k]), permutations = Nperm1)
        anosimSumm[k,u]=AA$delta
      }else{
        tmp.frame=as.data.frame(GPS.dist.cutree.rnd[,k])
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
ObsVsRnd=AllSummaryMean[,1]
submatrix.tmp=AllSummaryMean[1:(length(distcuts)-1),2:length(distcuts)]
Rnd.tmp=diag(submatrix.tmp)#,AllSummaryMean[length(distcuts)-1,1])
ObsVsRnd=cbind(ObsVsRnd,Rnd.tmp)
submatrix.tmp=AllSummaryStd[1:(length(distcuts)-1),2:length(distcuts)]
Rnd.tmp=diag(submatrix.tmp)#,AllSummaryStd[length(distcuts)-1,1])
ObsVsRnd=cbind(distcuts[1:length(distcuts)-1],ObsVsRnd,Rnd.tmp)
colnames(ObsVsRnd)=c("log10(distance)",paste(setMethod,"Obs",sep=" "),paste(setMethod,"Rnd Mean",sep=" "),"StdErr")

# Print output -----------------------------------------
#dir=dirOut
#setwd(dir)
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
# colnames(AllSummary)=c("distance",distcuts[1:length(distcuts)-1]) # and names

# ... Create file
Obs=AllSummaryMean[,1] # observations are in the first colum
AllSummaryTmp=AllSummaryMean # we want them instead in the diagonal
AllSummaryTmp[,1]=NA # remove them
diag(AllSummaryTmp)=Obs # and add them in the diag
AllSummaryTmp=t(AllSummaryTmp) # and transpose
AllSummary=as.data.frame(cbind(distcuts,AllSummaryTmp)) # add distances 
colnames(AllSummary)=c("distance",distcuts[1:length(distcuts)-1]) # and names


# ... prepare plot and go

shapes=seq(1:length(distcuts)-1)

df <- melt(AllSummary ,  id.vars = 'distance', variable.name = 'series')
figureOut=paste("DistanceDecay-Shuffled_",labelOut,".pdf",sep="")
colourSeries=c("black","blue","red","green4","grey",
               "cyan","orange","brown","violet","gold")
pdf(file=figureOut)
print(ggplot(df, aes(distance,value))+
  ylab(ylabel) +xlab("log10(distance) (meters)")+
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

# ... Store also this structure

write.table(AllSummary,sep="\t",file=fileMeanOut2,quote=FALSE,row.names=FALSE)

if(plotCtrl==1){
  ### Perform some more exploratory analysis
  mantel.pears=mantel(DistCor,log10(GPS.dist+1))
  mantel.spear=mantel(DistCor,log10(GPS.dist+1),method="spearman")
  
  # Retrieve results and plot, first pearson
  perm <- permustats(mantel.pears)
  outPlot=paste("mantelPears_logGPSdistVs",labelOut,".pdf",sep="")
  pdf(outPlot,height=7,width=14)
  g1 = densityplot(perm)
  g2 = qqmath(perm)
  grid.arrange(g1,g2,nrow = 1,top="Mantel test (Pearson)")
  dev.off()
  outFit=paste("Summary_mantelPears_logGPSdistVs",labelOut,".dat",sep="") # and the fit
  sink(outFit) # divert the output to the file
  print(mantel.pears) #  skip the limitation in lines s
  sink() # close the pipe
  
  # Retrieve results and plot, second Spearman
  perm <- permustats(mantel.spear)
  outPlot=paste("mantelSpearman_logGPSdistVs",labelOut,".pdf",sep="")
  pdf(outPlot,height=7,width=14)
  g1 = densityplot(perm)
  g2 = qqmath(perm)
  grid.arrange(g1,g2,nrow = 1,top="Mantel test (Spearman)")
  dev.off()
  outFit=paste("Summary_mantelSpearman_logGPSdistVs",labelOut,".dat",sep="") # and the fit
  sink(outFit) # divert the output to the file
  print(mantel.spear)
  sink() # close the pipe
}
