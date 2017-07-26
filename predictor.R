predictor=function(mirtron_input1,Models_all,windows){
  library("randomForest")
  library(MASS)
  library(ROCR)
  library(gridExtra)
  source("joyofdata.R")
  source('DataAnalysis_BK.R')
#load("\Files\Models.rda")
RFModels=list(Models_all[[1]][[3]],Models_all[[3]][[3]],Models_all[[5]][[3]],
              Models_all[[7]][[3]],Models_all[[9]][[3]],Models_all[[11]][[3]],
              Models_all[[13]][[3]])
probe_length=170
#itnumber=3
Results_all=list()
# dataset[[1]] are 5' tailed, dataset[[2]] are 3' tailed
datasets=list(intron_testing[which(intron_testing$class=="5p"),])#,intron_testing[which(intron_testing$class=="3p"),])
cat("Starting predictions...\n")
for (j in 1:length(windows)){
  cat(j,"...\n")
  window=windows[j]
  intron_testing=datasets[[1]]
  for (i in 1:dim(intron_testing)[1]){
    class=intron_testing$class[i]
    #gets apropriate slice of an intron
    if (class=="3p"){
      intron_testing$intronprobe[i]=substr(intron_testing$intron_seq[i],1,probe_length)
      
    } else if (class=="5p"){
      intron_testing$intronprobe[i]=substr(intron_testing$intron_seq[i],nchar(intron_testing$intron_seq[i])-probe_length+1,nchar(intron_testing$intron_seq[i]))
    } else {
      intron_testing$intronprobe[i]=intron_testing$intron_seq[i]
    }
  }
    x=intron_testing
    Scores=matrix(,ncol=1,nrow=dim(x)[1])
    types_num=Scores
  #generates subset (window) of intron probe sequence for every intron
  for (k in 1:dim(x)[1]){ #for every intron
    intronseq=data.frame(matrix(,nrow=nchar(x$intronprobe[k])-window+1,ncol=3))
    colnames(intronseq)=c("name","mature","position")
    
    for (i in 0:(nchar(x$intronprobe[k])-window)){ #for every subset
      intronseq$name[i+1]=x$hairpin_name[k]
      intronseq$mature[i+1]=substring(x$intronprobe[k],i+1,window+i) #window sequence
      #count position from the end of the intron
      if (x$class[k]=="3p" || x$class[k]=="mirtron"){ #if 3' tailed 
        intronseq$position[i+1]=i+1
      } else if (x$class[k]=="5p"){ #if 5' tailed
        intronseq$position[i+1]=nchar(x$intronprobe[k])-window-i+1
      } 
      
    }
    intronseq$mature_A=str_count(intronseq$mature,"A")*100/nchar(intronseq$mature)
    intronseq$mature_C=str_count(intronseq$mature,"C")*100/nchar(intronseq$mature)
    intronseq$mature_G=str_count(intronseq$mature,"G")*100/nchar(intronseq$mature)
    intronseq$mature_U=str_count(intronseq$mature,"U")*100/nchar(intronseq$mature)
    
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
    intronseq=intronseq[,-c(1,2)]
    #intronseq=data.frame(scale(intronseq))
    
    type=predict(RFModels[[j]],intronseq)
    types_num[k]=paste(type,collapse="")
    
    score=sum(as.integer(type)-1)
    Scores[k,1]=score
    #intronseq=rbind(intronseq,intronseq1) #final dataframe containing all windows
  }
    scores_number=6
    hist(Scores)
    
    Results=data.frame(matrix(,scores_number,2))
    colnames(Results)=c("Sens", "Spec")
  for (i in c(1:scores_number)){
      x$score=Scores
      x$answer=0
      x$answer[x$score>=i]=1
      x$answer[x$score<i]=0
      
      answer=table(x$answer,x$true)
      z=answer
      rfmSens=z[4]/(z[2]+z[4])
      rfmSpec=z[1]/(z[1]+z[3])
      Results$Sens[i]=rfmSens
      Results$Spec[i]=rfmSpec
      #cat("Senitivity equals ",rfmSens,"\n")
      #cat("Specificity equals ",rfmSpec,"\n")
  }
    x$type_num=as.character(types_num)
Results_all[[j]]=Results
}

lists=list(Results_all, x)
#print(Results)
return(lists)

}
