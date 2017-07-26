library("S4Vectors")
library("stats4")
library("IRanges")
library("XVector")
library("Biostrings")
library("stringr")
library("RRNA")
library("ggplot2")
setwd("/home/ror/Desktop/Mirtrony")
rm(list=ls())
source("overhangcount.R")
source("mirna_features.R")
source("intron_features.R")
source('DataAnalysis_BK.R')


startn=1
endn=224 # change it to find class=0 or 1
class="5p" #tailed
arm="3p"
window=25
mirratio=0.6
hairpin_length=120

mirtron_data=read.csv("./Data/mirtrondata.csv", header=TRUE, stringsAsFactors = FALSE)
mirna_input1=data.frame(hairpin_seq=mirtron_data$hairpin_seq, db=mirtron_data$dotbracket, fe=mirtron_data$fe, 
                        mature5p_seq=mirtron_data$mature5p_seq, mature3p_seq=mirtron_data$mature3p_seq, stringsAsFactors = FALSE)
mirna_5p=mirna_input1[c(which(mirtron_data$class=="5p")),]
mirna_3p=mirna_input1[c(which(mirtron_data$class=="3p")),]
mirna_mirtron=mirna_input1[c(which(mirtron_data$class=="mirtron")),]
#mirtron_mirna=mirna_features(mirna_input1, random=FALSE) #cal, mirna features
#mirtron_mirna$class=1 #1 is for mirtron
#mirtron_mirna$type=mirtron_data$class

#introns=introns(mirtron_data)
intron_data=read.csv("./Data/introns.csv", header=TRUE, stringsAsFactors = FALSE)
intron_input1=data.frame(hairpin_seq=intron_data$intron_seq, 
                        intron5p_seq=intron_data$intron_seq, intron3p_seq=intron_data$intron_seq, stringsAsFactors = FALSE)
intron_5p=intron_input1[c(which(mirtron_data$class=="5p")),]
intron_3p=intron_input1[c(which(mirtron_data$class=="3p")),]
intron_mirtron=intron_input1[c(which(mirtron_data$class=="mirtron")),]

#source("otherintrons.R")
intron_data2=read.csv("./Data/all_other_introns.csv", header=TRUE, stringsAsFactors = FALSE)
intron_input2=data.frame(hairpin_seq=intron_data2$intron_seq, 
                         intron5p_seq=intron_data2$intron_seq, intron3p_seq=intron_data2$intron_seq, stringsAsFactors = FALSE)
#intron_intron$class=0

#intron=intron_input1[startn:endn,]
#mirtron_data_i=mirtron_data[startn:endn,]
#mirtron_mirna_i=mirtron_mirna[startn:endn,]
intron$intron5p_seq=chartr("T","U",intron$intron5p_seq)
intron$intron3p_seq=chartr("T","U",intron$intron3p_seq)


if (class=="3p"){
  hairpin=data.frame(hairpin=intron$intron5p_seq[which(mirtron_data_i$class=="3p")])
  mirtron1=data.frame(mirtron_data_i[which(mirtron_data_i$class=="3p"),])
  mirtron2=data.frame(mirtron_mirna_i[which(mirtron_data_i$class=="3p"),])
  
  } else if (class=="5p"){
    hairpin=data.frame(hairpin=intron$intron3p_seq[which(mirtron_data_i$class=="5p")])
    mirtron1=data.frame(mirtron_data_i[which(mirtron_data_i$class=="5p"),])
    mirtron2=data.frame(mirtron_mirna_i[which(mirtron_data_i$class=="5p"),])
    
  } else{
    hairpin=data.frame(hairpin=intron$intron3p_seq[which(mirtron_data_i$class=="mirtron")])
    mirtron1=data.frame(mirtron_data_i[which(mirtron_data_i$class=="mirtron"),])
    mirtron2=data.frame(mirtron_mirna_i[which(mirtron_data_i$class=="mirtron"),])
}
if (arm=="3p"){
  x=substr(hairpin$hairpin,ceiling(nchar(as.character(hairpin$hairpin))/3)+1,nchar(as.character(hairpin$hairpin)))
  mature=data.frame(mature=mirtron1$mature3p_seq)
  position=data.frame(position=mirtron2$mature3pposition)
} else{
  x=substr(hairpin$hairpin,1,ceiling(nchar(as.character(hairpin$hairpin))*2/3))
  mature=data.frame(mature=mirtron1$mature5p_seq)
  position=data.frame(position=mirtron2$mature5pposition)
}

intronseq=matrix(,nrow=1,ncol=3)
colnames(intronseq)=c("mature","position","type")
known=data.frame(mature=mature,position=position,class=class,arm=arm,type=1)

for (k in 1:dim(hairpin)[1]){
  intronseq1=data.frame(matrix(,nrow=nchar(x[k])-window,ncol=3))
  colnames(intronseq1)=c("mature","position","type")
  for (i in 0:(nchar(x[k])-window)){
    intronseq1[i+1,1]=substring(x[k],i+1,window+i)
    intronseq1[i+1,3]=0
    if (arm=="3p"){
      intronseq1[i+1,2]=nchar(x[k])-window-i
    } else {
      intronseq1[i+1,2]=i+1
    }
   
  }
  intronseq=rbind(intronseq,intronseq1)
  
}
known=known[,c(-3,-4)]
known$mature=as.character(known$mature)
intronseq=intronseq[-1,]
intronseq$mature=as.character(intronseq$mature)
C=c()
for (i in 1:dim(known)[1]){
  z=known$mature[i]
  mirwindow=ceiling(nchar(z)*mirratio)
  k=1
  while (k+mirwindow<=nchar(z)){
    c=which(!is.na(str_locate(intronseq$mature,substring(z,k,k+mirwindow)))[,1])
    C=c(c,C)
    k=k+1
  }
  
}
if (!is.null(C)){
  #intronseq=intronseq[-C,]
}

  

input=rbind(intronseq,known)
input$mature=as.character(input$mature)
results=data.frame(position=as.numeric(input$position))
results$mature_length=nchar(input$mature)
results$mature_A=str_count(input$mature,"A")*100/nchar(input$mature)
results$mature_C=str_count(input$mature,"C")*100/nchar(input$mature)
results$mature_G=str_count(input$mature,"G")*100/nchar(input$mature)
results$mature_U=str_count(input$mature,"U")*100/nchar(input$mature)
results$seq=input$mature
results$class=input$type
results$class[C]=1########

pcdata5=results
pcdata5=pcdata5[,c(-2,-7,-8)]
pca5=prcomp(pcdata5, retx=TRUE, center=TRUE, scale=TRUE)
row.names(pca5$x)=1:dim(pcdata5)[1]

labels=factor(c(rep('introns',dim(intronseq)[1]),rep('mirtrons',dim(pcdata5)[1]-dim(intronseq)[1])),levels=c('introns','mirtrons'))
gg2_pca5<-BK_gg2_biplot(pca5,gr.shape = labels, gr.colour =labels, size=2, attr_text_size=6,
                       limits = list(x=c(-5,5),y=c(-5,5)))

title=paste(paste(paste(class, "tailed",sep="-"),arm,sep=", "), "arm")
  
plot(gg2_pca5, main=title)

