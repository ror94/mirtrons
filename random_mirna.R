library("Biostrings")
random_mirna=data.frame(matrix(,nrow=200,ncol=2))
colnames(random_mirna)=c("mature5p_seq","mature3p_seq")

for (i in 1:dim(random_mirna)[1]) {
  length5p=sample(17:27,1) #losowa długość miRNA
  length3p=sample(17:27,1)
  seq5p=""
  seq3p=""
  for (k in 1:length5p){
    base5p=sample(c("A","C","G","U"),1) #losowa sekwencja 5p
    seq5p=paste(seq5p,base5p,sep="")
  }
  overhang=sample(-5:5,1) #losowy overhang (teoretycznie nie powinien być rozkład jednostajny)
  overhangseq=""          # dodatni overhang tzn ze wystaje ramie 3p, ujemny 5p

  
  if (overhang<0){ 
    for (l in -1:overhang){
      baseoh=sample(c("A","C","G","U"),1) 
      overhangseq=paste(overhangseq,baseoh,sep="") #losowa sekwencja overhangu
    }
  # tworze sekwencję komplementarną z przesunięciem
  seq3p=paste(paste(reverseComplement(RNAString(substring(seq5p,1,(nchar(seq5p)+overhang))))),overhangseq,sep="")
  } 
  else if (overhang==0) {
    seq3p=paste(reverseComplement(RNAString(seq5p))) #tworze sekwencję komplementarną
  } 
  else {
    for (l in 1:overhang){
      baseoh=sample(c("A","C","G","U"),1)
      overhangseq=paste(overhangseq,baseoh,sep="")
  
    }
    #tworzę sekwencję komplementarną z przesunięciem w drugą stronę
    seq3p=paste(reverseComplement(RNAString(seq5p)))
    seq5p=paste(overhangseq,substring(seq5p, 1, nchar(seq5p)-overhang),sep="")
  }
  
 
 random_mirna$mature5p_seq[i]=seq5p
 random_mirna$mature3p_seq[i]=seq3p
 random_mirna$oh[i]=overhang
}

write.csv(random_mirna, "./Data/random_mirna.csv")
