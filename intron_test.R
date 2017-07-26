if (Sys.info()[[1]]=="Windows"){
  setwd("C:/Users/Ror/Desktop/R Mirtrony/Mirtrony")
} else{
  setwd("/home/ror/Desktop/R mirtrony/Mirtrony")
}
rm(list=ls())
library(e1071)
library("S4Vectors")
library("stats4")
library("IRanges")
library("XVector")
library("Biostrings")
library("stringr")
library("RRNA")
library("ggplot2")
library("rpart")
library("randomForest")
library(AUC)
library(ROCR)
library(rpart)
library(psych)
library(tree)
library(MASS)
source("overhangcount.R")
source("mirna_features.R")
source("intron_features.R")
source('DataAnalysis_BK.R')
source('LogReg.R')
source('predictor.R')
source("BinEval.R")



#MIRTRONS
mirtron_data=read.csv("./Data/mirtrondata.csv", header=TRUE, stringsAsFactors = FALSE)
mirna_input1=data.frame(hairpin_name=mirtron_data$hairpin_name,hairpin_seq=mirtron_data$hairpin_seq, db=mirtron_data$dotbracket, fe=mirtron_data$fe, 
                        mature5p_seq=mirtron_data$mature5p_seq, mature3p_seq=mirtron_data$mature3p_seq, stringsAsFactors = FALSE)


#INTRONS CODING MIRNA
intron_data=read.csv("./Data/introns.csv", header=TRUE, stringsAsFactors = FALSE)
intron_input1=data.frame(intron_seq=intron_data$intron_seq, class=mirtron_data$class, stringsAsFactors = FALSE)
intron_input1$intron_seq=chartr("T","U",intron_input1$intron_seq)


#INTRONS NOT CODING MIRNA
intron_data2=read.csv("./Data/negativeintrons.csv", header=TRUE, stringsAsFactors = FALSE)
names(intron_data2)[1]='hairpin_name'
intron_data2$hairpin_seq=""
intron_data2$db=""
intron_data2$fe=""
intron_data2$mature5p_seq=""
intron_data2$mature3p_seq=""
intron_data2$class="5p"
intron_data2$intron_seq=chartr("T","U",intron_data2$intron_seq)


intron_data3=intron_data2
intron_data3$class="3p"
intron_data2=rbind(intron_data2,intron_data3)
intron_data2$true=0

#adding intron sequences to miRNA data
mirna_input1=cbind(mirna_input1,intron_input1)
mirna_input1$true=1
mirna_input1=rbind(mirna_input1,intron_data2) #ready list
mirna_input2=mirna_input1#[which(mirna_input1$class=="5p"),]
X=BK_sample_crossVal(dim(mirna_input2)[1],3)
intron_testing=mirna_input2[X[[1]],]
mirna_input1=mirna_input2[-X[[1]],]
mirna_input1$true=NULL


#windows=seq(7,31,4)
windows=c(23)
probe_length=170

Results_all=list()
Models_all=list()

for (j in 1:length(windows)){
    #
    window=windows[j]  
    for (i in 1:dim(mirna_input1)[1]){
      class=mirna_input1$class[i]
      #gets apropriate slice of an intron
      if (class=="3p"){
        mirna_input1$intronprobe[i]=substr(mirna_input1$intron_seq[i],1,probe_length)
        
      } else if (class=="5p"){
        mirna_input1$intronprobe[i]=substr(mirna_input1$intron_seq[i],nchar(mirna_input1$intron_seq[i])-probe_length+1,nchar(mirna_input1$intron_seq[i]))
      } else {
        mirna_input1$intronprobe[i]=mirna_input1$intron_seq[i]
      }
      #locates 5p and 3p arm
      loc5p=str_locate(mirna_input1$intronprobe[i],mirna_input1$mature5p_seq[i])
      loc3p=str_locate(mirna_input1$intronprobe[i],mirna_input1$mature3p_seq[i])
      x=as.character(rep(0,nchar(mirna_input1$intronprobe[i])))
      #puts 1 where 5p arm, 2 where 3p arm nad 3 where in between (loop)
      if (loc5p[2]-loc5p[1] != 0 & loc3p[2]-loc3p[1] != 0){
      x=replace(x,loc5p[1]:loc5p[2],1)
      x=replace(x,loc3p[1]:loc3p[2],2)
      x=replace(x,(loc5p[2]+1):(loc3p[1]-1),3)
      }
      x=paste(x,collapse="")
      mirna_input1$numerical[i]=x
      
    }
    #merges data of introns containing miRNAs (mirna_input1) 
    #and introns not containing miRNAs (intron_data2)
    x=mirna_input1
    intronseq=matrix(,nrow=0,ncol=5)
    colnames(intronseq)=c("name", "mature","position","type","class")
    
    #generates subset (window) of intron probe sequence for every intron
    for (k in 1:dim(x)[1]){ #for every intron
      intronseq1=data.frame(matrix(,nrow=nchar(x$intronprobe[k])-window+1,ncol=5))
      colnames(intronseq1)=c("name","mature","position","type","class")
      for (i in 0:(nchar(x$intronprobe[k])-window)){ #for every subset
        intronseq1$name[i+1]=x$hairpin_name[k]
        intronseq1$mature[i+1]=substring(x$intronprobe[k],i+1,window+i) #window sequence
        #type of nucleotide in the middle (0,1,2 or 3)
        intronseq1$type[i+1]=substr(x$numerical[k],i+ceiling(window/2),i+ceiling(window/2))
        #count position from the end of the intron
        if (x$class[k]=="3p" || x$class[k]=="mirtron"){ #if 3' tailed 
          intronseq1$position[i+1]=i+1
        } else if (x$class[k]=="5p"){ #if 5' tailed
          intronseq1$position[i+1]=nchar(x$intronprobe[k])-window-i+1
        } 
        intronseq1$class=x$class[k]
        
      }
      intronseq=rbind(intronseq,intronseq1) #final dataframe containing all windows
     }
    names=intronseq[,1] 
    intronseq=intronseq[,-1] #delet names column
    #translates numerical values to  character
    intronseq$type=replace(intronseq$type,which(intronseq$type=="0"),"intron")
    intronseq$type=replace(intronseq$type,which(intronseq$type=="1"),"5p arm")
    intronseq$type=replace(intronseq$type,which(intronseq$type=="2"),"3p arm")
    intronseq$type=replace(intronseq$type,which(intronseq$type=="3"),"loop")
    #calculates precentage
    intronseq$mature_A=str_count(intronseq$mature,"A")*100/nchar(intronseq$mature)
    intronseq$mature_C=str_count(intronseq$mature,"C")*100/nchar(intronseq$mature)
    intronseq$mature_G=str_count(intronseq$mature,"G")*100/nchar(intronseq$mature)
    #intronseq$mature_U=str_count(intronseq$mature,"U")*100/nchar(intronseq$mature)
    
    intronseq$mature_AA=str_count(intronseq$mature,"AA")
    intronseq$mature_AU=str_count(intronseq$mature,"AU")
    intronseq$mature_AC=str_count(intronseq$mature,"AC")
    intronseq$mature_AG=str_count(intronseq$mature,"AG")
    
    intronseq$mature_UA=str_count(intronseq$mature,"UA")
    intronseq$mature_UU=str_count(intronseq$mature,"UU")
    intronseq$mature_UC=str_count(intronseq$mature,"UC")
    intronseq$mature_UG=str_count(intronseq$mature,"UG")
    
    intronseq$mature_CA=str_count(intronseq$mature,"CA")
    intronseq$mature_CU=str_count(intronseq$mature,"CU")
    intronseq$mature_CC=str_count(intronseq$mature,"CC")
    intronseq$mature_CG=str_count(intronseq$mature,"CG")
    
    intronseq$mature_GA=str_count(intronseq$mature,"GA")
    intronseq$mature_GU=str_count(intronseq$mature,"GU")
    intronseq$mature_GC=str_count(intronseq$mature,"GC")
    intronseq$mature_GG=str_count(intronseq$mature,"GG")
    
    #sequences of 5'-tailed introns
    seq5prim=intronseq[which(intronseq$class=="5p"),]
    type5prim=seq5prim$type
    #sequence of 5'tailed introns (without windows classified as introns)
    #seq5prim_intron=seq5prim[-c(which(type5prim=="intron")),]
    #type5prim_intron=seq5prim_intron$type
    #sequences of 3' tailed introns
    seq3prim=intronseq[which(intronseq$class=="3p"),]
    type3prim=seq3prim$type
    #sequences od 3'failed intron (without windows classified as introns)
    #seq3prim_intron=seq3prim[-c(which(type3prim=="intron")),]
    #type3prim_intron=seq3prim_intron$type
    #sequences of mirtrons
    seqmirtron=intronseq[which(intronseq$class=="mirtron"),]
    typemirtron=seqmirtron$type
    
    seq5prim=seq5prim[,-c(1,3,4)]
    seq3prim=seq3prim[,-c(1,3,4)]
    seqmirtron=seqmirtron[,-c(1,3,4)]
    
    
    pca=prcomp(seq5prim, retx=TRUE, center=TRUE, scale=TRUE)
    labels=factor(type5prim,levels=c('5p arm','3p arm','loop','intron'))
    gg2_pca<-BK_gg2_biplot(pca,gr.shape = labels, gr.colour =labels, size=2, attr_text_size=4,
                            limits = list(x=c(-5,6),y=c(-5,5)))
    plot(gg2_pca)
    
    pca2=prcomp(seq3prim, retx=TRUE, center=TRUE, scale=TRUE)
    labels=factor(type3prim,levels=c('5p arm','3p arm','loop','intron'))
    gg2_pca2<-BK_gg2_biplot(pca2,gr.shape = labels, gr.colour =labels, size=2, attr_text_size=4,
                           limits = list(x=c(-5,6),y=c(-5,5)))
    #plot(gg2_pca2)
    
    
    pca3=prcomp(seqmirtron, retx=TRUE, center=TRUE, scale=TRUE)
    labels=factor(typemirtron,levels=c('5p arm','3p arm','loop','intron'))
    gg2_pca3<-BK_gg2_biplot(pca3,gr.shape = labels, gr.colour =labels, size=2, attr_text_size=4,
                           limits = list(x=c(-5,6),y=c(-5,5)))
    #plot(gg2_pca3)
    
    
    
    ####################################################################################################################################
    #pcas=list(pca,pca2)
    
    #pcas=list(seq5prim,seq3prim)
    pcas=list(seq5prim)
    #types=list(type5prim,type3prim)
    types=list(type5prim)
    rm(intronseq,gg2_pca,gg2_pca2,gg2_pca3,seq5prim,seq3prim,seqmirtron,type5prim,type3prim,typemirtron,pca,pca2,pca3)
    tailed=c('5\'-tailed','3\'-tailed')
    
    for (k in 1:length(pcas)){
      #Z=data.frame(pcas[[k]]$x[,1:15])
      Z=data.frame(pcas[[k]])
      #means=apply(Z[1:21],2,mean)
      #stds=apply(Z[1:21],2,sd)
      #Z=data.frame(scale(Z))
      Z$class=types[[k]]
      
      Z_intron=Z
      Z_intron$class[which(Z$class=="intron")]=1
      Z_intron$class[which(Z$class!="intron")]=0
      Z_intron$class=as.numeric(Z_intron$class)
      
      Z_5parm=Z
      Z_5parm$class[which(Z$class=="5p arm")]=1
      Z_5parm$class[which(Z$class!="5p arm")]=0
        Z_5parm$class=as.numeric(Z_5parm$class)
      
      Z_3parm=Z
      Z_3parm$class[which(Z$class=="3p arm")]=1
      Z_3parm$class[which(Z$class!="3p arm")]=0
      Z_3parm$class=as.numeric(Z_3parm$class)
      
      Z_loop=Z
      Z_loop$class[which(Z$class=="loop")]=1
      Z_loop$class[which(Z$class!="loop")]=0
      Z_loop$class=as.numeric(Z_loop$class)
      Datasets=list(Z_intron, Z_5parm,Z_3parm,Z_loop)
      rm(Z,Z_5parm,Z_3parm,Z_intron,Z_loop)
      
      Results=list()
      Models=list()
      #colnames(Results)=c('LRSens','LRSpec','RFSens','RFSpec',"LDASens",'LDASpec')
      #rownames(Results)=c('Intron','5p arm', '3p arm', 'loop')
      for(i in 1:length(Datasets)) {
        cat(c('Learning',rownames(Results)[i],'of',tailed[k],'mirtron with window =',windows[j],'...\n'))
        list1=LogReg(Datasets[[i]],3)
        print(list1[[1]])
        Results[[i]]=list1[[1]]
        Models[[i]]=list1[[2]]
      }
      Results_all[[(j-1)*2+k]]=rbind(Results[[1]],Results[[2]],Results[[3]],Results[[4]])
      Models_all[[(j-1)*2+k]]=Models
      names(Results_all)[(j-1)*2+k]=paste(paste('window=',windows[j]),paste(tailed[k]) )
      names(Models_all)[(j-1)*2+k]=paste(paste('window=',windows[j]),paste(tailed[k]) )
       }
    
}
lapply(Results_all, function(x) write.table( data.frame(x), 'test.csv'  , append= T, sep=',' ))
RFModels=list(Models_all[[11]][[3]])
cat("Saving Models...\n")
save(RFModels, file="./Files/RFModels.rda")
save(Models_all, file="./Files/Models.rda")
cat("Running predictior...\n")
#final=predictor(intron_testing,Models_all,windows=15)
#write.csv(final[[1]],"./Files/final.csv")
#write.csv(final[[2]],"./Files/final2.csv")
#cat("Done")
                                                                        

