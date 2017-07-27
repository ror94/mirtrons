rm(list=ls())
if (Sys.info()[[1]]=="Windows"){
  setwd("D:/BACKUP/Workspace/R/Mirtrons")
} else{
  setwd("~/Workspace/R/Mirtrons")
}
library(Boruta) #rf importance
library(agricolae) #tukey test
library("ggplot2")
library(psych) # pair_panels
library(ggbiplot)
library("stringr") #str_find
library("RRNA") #makect
library(FactoMineR) #PCA
library("factoextra")
library("randomForest")
library(ROCR) #roc curve
library(MASS) #lda
library(tree) #tree
library(e1071) #svm
#library(party) #cforest
library(caret) # RFE
library(mlr)
library(TunePareto) # stratified

source("loopscount.R")
source("mirna_features.R")
source("overhangcount.R")
source("LogReg.R")
source("DataAnalysis_BK.R")
source("BinEval.R")

#library("rpart")
#library(AUC)
#library("S4Vectors")
#library("stats4")
#library("IRanges")
#library("XVector")
#library("Biostrings")

#####################################################################################################################################################################

#CANONICAL MIRNAs
canonical_data=read.csv("./Data/prep_names.csv", header=TRUE, stringsAsFactors = FALSE)
canonical_data=canonical_data[-c(380,702,720,813,889),]
mirna_input2=data.frame(hairpin_seq=canonical_data$hairpin_seq, db=canonical_data$dotbracket, fe=canonical_data$fe, 
                        mature5p_seq=canonical_data$mature5p_seq, mature3p_seq=canonical_data$mature3p_seq, stringsAsFactors = FALSE)
canonical_mirna=mirna_features(mirna_input2,random=FALSE)
canonical_mirna$class=0
#MIRTRONS
mirtron_names=read.table("./Data/mirtron_names.txt")
Index=c()
for (i in 1:dim(mirtron_names)[1]){
    index=grep(paste("^",mirtron_names[i,1],"$",sep=""),canonical_data$hairpin_name)
    Index=rbind(Index,index)
}
mirtron_mirna=canonical_mirna[Index,]
canonical_mirna=canonical_mirna[-Index,]
mirtron_mirna$class=1


#TEST
test_data=read.csv("./Data/testdata.csv", header=TRUE, stringsAsFactors = FALSE)
test_data=test_data[-c(1,22,103,139,151,164,165,182,202),]
test_input1=data.frame(hairpin_seq=test_data$hairpin_seq, db=test_data$dotbracket, fe=test_data$fe, 
                        mature5p_seq=test_data$mature5p_seq, mature3p_seq=test_data$mature3p_seq, stringsAsFactors = FALSE)
test_mirna=mirna_features(test_input1,random=FALSE)
test_mirna$class=1 #theoretically mirtrons


#RANDOM
random_data=read.csv("./Data/randomdata.csv", header=TRUE, stringsAsFactors = FALSE)
random_input=data.frame(hairpin_seq="T", mature5p_seq=random_data$mature5p_seq, mature3p_seq=random_data$mature3p_seq, stringsAsFactors = FALSE)

random_mirna=mirna_features(random_input,random=TRUE)
random_mirna$class=1


######################################################################################################
#Plots
source('mirnaplots.R')
#mirtronplots_mirna=mirnaplots(mirtron_mirna)
#canonicalplots_mirna=mirnaplots(canonical_mirna)
#testplots_mirna=mirnaplots(test_mirna)

P.values=data.frame(matrix(,nrow=dim(canonical_mirna)[2]-1,ncol=5))
sign=c('Not Significant','Significant')
for (i in 1:(ncol(canonical_mirna)-1)){
  P.values[i,1]=wilcox.test(mirtron_mirna[,i],canonical_mirna[,i],alternative="two.sided")$p.value
  P.values[i,2]=ks.test(mirtron_mirna[,i],canonical_mirna[,i],alternative="two.sided")$p.value
  
  tukey=data.frame(input=c(mirtron_mirna[,i],canonical_mirna[,i]),classes=c(rep(1,dim(mirtron_mirna)[1]),rep(0,dim(canonical_mirna)[1])))
  means=tapply(tukey$input,tukey$classes, mean)
  ajuste <- lm( tukey$input ~ tukey$classes)
  tukey.alpha=0.005
  h=HSD.test(ajuste, 'tukey$classes',alpha=tukey.alpha)
  P.values[i,3]=sign[length(levels(h$groups$M))]
  P.values[i,4]=means[1]
  P.values[i,5]=means[2]
}
rownames(P.values)=colnames(canonical_mirna[,1:(dim(canonical_mirna)[2]-1)])
colnames(P.values)=c('Wilcoxon','Kolmogorov-Smirnov',paste('Tukey',tukey.alpha),'Mirtron mean','Canonical mean')
cat("\nStatistical tests\n")
print(P.values)

####################################################################################################################################################################
#PCA
canonical_mirna$mature3pposition=NULL
canonical_mirna$mature5pposition=NULL
canonical_mirna$mature5p_U=NULL
canonical_mirna$mature3p_U=NULL
canonical_mirna$hairpin_U=NULL

mirtron_mirna$mature3pposition=NULL
mirtron_mirna$mature5pposition=NULL
mirtron_mirna$mature5p_U=NULL
mirtron_mirna$mature3p_U=NULL
mirtron_mirna$hairpin_U=NULL

random_mirna$mature5p_U=NULL
random_mirna$mature3p_U=NULL
random_mirna$hairpin_U=NULL

pcdata_ml=rbind(mirtron_mirna,canonical_mirna)

#pairs.panels(pcdata_ml)
# PCA canonical vs mirtron
pcdata=pcdata_ml[,-ncol(pcdata_ml)]
pca=prcomp(pcdata, retx=TRUE, center=TRUE, scale=TRUE)
labels=factor(c(rep('mirtron',dim(mirtron_mirna)[1]),rep('canonical',dim(pcdata)[1]-dim(mirtron_mirna)[1])),levels=c('mirtron','canonical'))
g <- ggbiplot(pca, obs.scale = 1, var.scale = 1, 
              groups = labels, ellipse = F, 
              circle = F)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g) #mirtron vs canonical

test_mirna$mature3pposition=NULL
test_mirna$mature5pposition=NULL
test_mirna$mature5p_U=NULL
test_mirna$mature3p_U=NULL
test_mirna$hairpin_U=NULL


# Added test
pcdata5=rbind(pcdata,test_mirna[,-ncol(test_mirna)])
pca5=pca
pca5$x=scale(pcdata5, pca$center, pca$scale) %*% pca$rotation
rownames(pca5$x)=1:dim(pcdata5)[1]
labels=factor(c(rep('mirtron',dim(mirtron_mirna)[1]),rep('canonical',dim(pcdata5)[1]-dim(test_mirna)[1]-dim(mirtron_mirna)[1]),
                rep('test',dim(pcdata5)[1]-dim(canonical_mirna)[1]-dim(mirtron_mirna)[1]))
              ,levels=c('mirtron','canonical','test'))
g2 <- ggbiplot(pca5, obs.scale = 1, var.scale = 1, 
              groups = labels, ellipse = F, 
              circle = F)

g2 <- g2 + scale_color_discrete(name = '')
g2 <- g2 + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g2)



#pca_new=PCA(pcdata)
#f1=fviz_contrib(pca_new,choice = "var", axes = 1)
#f2=fviz_contrib(pca_new,choice = "var", axes = 2)
#f1_2=fviz_contrib(pca_new,choice = "var", axes = 1:2)
#plot(f1)
#plot(f2)
#plot(f1_2)
#fviz_pca_var(pca_new, col.var="contrib") +
#  scale_color_gradient2(low="black", mid="blue", 
#                        high="red", midpoint=6) + theme_minimal)

# 3d PCA
library(rgl)
colors = replace(pcdata_ml$class, with(pcdata_ml, which(class == 1)), "blue")
colors = replace(colors,which(colors== 0), "green")
coords <- NULL
for (i in 1:nrow(pca$rotation)) {
  coords <- rbind(coords, rbind(c(0,0,0),pca$rotation[i,1:3]))
}
plot3d(pca$x[,1:3], col=colors, size = 5)
text3d(pca$rotation[,1:3]*10, texts=rownames(pca$rotation), col="red")
lines3d(coords*10, col="red", lwd=2)

####################################################################################################
# Classification

source("LogReg.R")
cat("\n1x5-fold Cross-validation classification\n")
x=LogReg(pcdata_ml, itnumber=10)
print(x[[1]])
pred=predict(x[[3]],test_mirna) #predict on svm model
print(table(pred))
mean(as.double(pred)-1)

####################################################################################################
#Single feature

Singlef=data.frame()
for (i in 1:(dim(pcdata_ml)[2]-1)){
  singlef=LogReg(pcdata_ml[,c(i,dim(pcdata_ml)[2])],3)
  Singlef=rbind(Singlef,singlef[[1]][5,])
}
rownames(Singlef)=colnames(pcdata)
colnames(Singlef)=colnames(singlef[[1]])
ord=with(Singlef,order(-AUC))
Singlef=Singlef[ord,]
print(Singlef)

####################################################################################################
#Boruta & RFE
classes = replace(pcdata_ml$class, with(pcdata_ml, which(class == "1")), "mirtron")
classes = replace(classes, which(classes == "0"), "canonical")
svmProfile <-rfe(pcdata, as.factor(classes),sizes=c(1:(ncol(canonical_mirna)-1)),
                rfeControl = rfeControl(functions = caretFuncs,
                                        verbose = FALSE, method = "cv", number = 3),
                method = "svmRadial")
plot(svmProfile)
predictors(svmProfile)
plot(svmProfile, type=c("g", "o"))

imp2=getImpRfGini(pcdata_ml[,-(ncol(canonical_mirna))],pcdata_ml$class)
ord=order(-as.double(imp2))
Imp=data.frame(RFBoruta=names(imp2)[order(-as.double(imp2))][1:length(predictors(svmProfile))],RFrfe = predictors(svmProfile))
cat("\nImportance by random forest using Boruta package and RFE\n")
print(Imp)

####################################################################################################
#Classification using top features

top_feat = cbind(pcdata_ml[,predictors(svmProfile)],class = pcdata_ml$class)
cat("\n1x5-fold Cross-validation classification\n")
x=LogReg(top_feat, itnumber=10)
print(x[[1]])
pred=predict(x[[3]],test_mirna) #predict on svm model
print(table(pred))
mean(as.double(pred)-1)





