

mature <- function(data){
number=nrow(data)
mature = readRNAStringSet("./Fasta/mature.fa")
sequences=matrix(,number,2)
data$hairpin_seq=chartr("T","U",data$hairpin_seq)

for(i in 1:number){
  name=data$hairpin_name[i]
  name5p=paste(name,"-5p",sep="")
  name3p=paste(name,"-3p",sep="")
  index5p=which(!is.na(str_locate(names(mature),name5p)[,1]), TRUE)
  index3p=which(!is.na(str_locate(names(mature),name3p)[,1]), TRUE)
  index=which(!is.na(str_locate(names(mature), paste(name,""))[,1]), TRUE)
  
  if (!isEmpty(index3p)&!isEmpty(index5p)){
    sequences[i,1]=paste(mature)[index5p]
    sequences[i,2]=paste(mature)[index3p]
  } else if(!isEmpty(index5p)){
    sequences[i,1]=paste(mature)[index5p]
    sequences[i,2]=""
  } else if (!isEmpty(index3p)){
  sequences[i,1]=""
  sequences[i,2]=paste(mature)[index3p]
  } else if (!isEmpty(index)&isEmpty(index5p)&isEmpty(index3p)) {
    hairpin=data$hairpin_seq[i]
    loc=str_locate(hairpin,paste(mature)[index])
    if(loc[1]<nchar(hairpin)/2){
      sequences[i,1]=paste(mature)[index]
      sequences[i,2]=""
      }else {
      sequences[i,2]=paste(mature)[index]
      sequences[i,1]=""
    }

  }
  else{
    #sequences[i,1]=paste(mature)[index5p]
    #sequences[i,2]=paste(mature)[index3p]
  }
}
sequences=data.frame(sequences)
colnames(sequences)=c("mature5p_seq", "mature3p_seq")
sequences$mature5p_seq=as.character(sequences$mature5p_seq)
sequences$mature3p_seq=as.character(sequences$mature3p_seq)
sequences$mature5p_seq[is.na(sequences$mature5p_seq)]=""
sequences$mature3p_seq[is.na(sequences$mature3p_seq)]=""
return(sequences)
}
