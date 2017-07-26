library("S4Vectors")
library("stats4")
library("IRanges")
library("XVector")
library("Biostrings")
library("stringr")
setwd('C:/Users/Ror/Desktop/R mirtrony/Mirtrony/Prep_new')
precursors=data.frame(read.csv('hsa.csv', header=TRUE, sep=";"))
mature=read.csv('hsa2.csv', header=TRUE, sep=";")
high_conf=read.csv('high_conf.csv', header=TRUE, sep=";")
positive=read.csv('mirtrons.csv', header=TRUE, sep=";")
arms=matrix(,nrow=length(precursors$Chromosome),ncol=6)
for (i in 1:length(precursors$Chromosome)){
  indices=which(mature$Derives_from==precursors$ID[i])
  for (k in 1:length(indices)){
    if (mature$Arm[indices[k]]=="5p"){
      arms[i,1]=as.character(mature$Name[indices[k]])
    } else if (mature$Arm[indices[k]]=="3p") {
      arms[i,2]=as.character(mature$Name[indices[k]])
    } else {
      arms[i,3]=as.character(mature$Name[indices[k]])
    }
    if (length(indices)==2){
      arms[i,4]=1
    }
    if (precursors$ID[i] %in% high_conf$ID){
      arms[i,5]=1
    } else{
      arms[i,5]=0
    }
    if (precursors$Name[i] %in% positive$Name){
      arms[i,6]=as.character(positive$Type[which(as.character(precursors$Name[i])==as.character(positive$Name))])
    } else{
      arms[i,6]="canonical"
    }
  }
}

precursors$arm5p=arms[,1]
precursors$arm3p=arms[,2]
precursors$dk=arms[,3]
precursors$both=arms[,4]
precursors$high_conf=arms[,5]
precursors$mirtron=arms[,6]
precursors_both=precursors[which(precursors$both==1),-c(2,10,11)]
#precursors_conf=precursors_both[which(precursors_both$high_conf==1),-9]


library(Biostrings)
s = readRNAStringSet("mature.fa")
nam=names(s)
test_data=read.csv("../Data/testdata.csv", header=TRUE, stringsAsFactors = FALSE)
X=matrix(,ncol=3)
Y=matrix(,ncol=3)
for (i in 1: length(test_data$mature5p_seq)){
  seq5p=test_data$mature5p_seq[i]
  seq3p=test_data$mature3p_seq[i]
  x=which(!is.na(str_locate(as.character(s),seq5p)))
  y=which(!is.na(str_locate(as.character(s),seq3p)))
  lx=length(x)/2
  ly=length(y)/2
  if (lx){
  X=rbind(X,c(i,test_data$hairpin_name[i],names(s)[x[1]]))
  }
  if (ly){
  Y=rbind(Y,c(i,test_data$hairpin_name[i],names(s)[y[1]]))
  lx=0
  ly=0
  }
}
X=X[-1,]
Y=Y[-1,]

