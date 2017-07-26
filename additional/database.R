rm(list = ls())
setwd("C:/Users/Ror/Desktop/R Mirtrony/Mirtrony")

prep=read.csv("./Prep_new/prep.csv", header=TRUE, stringsAsFactors = FALSE)
s=readRNAStringSet("./Fasta/hairpin.fa")
m=readRNAStringSet("./Fasta/mature.fa")

mirnames=substr(names(s),1,regexpr(" ",names(s))-1)
maturenames=substr(names(m),1,regexpr(" ",names(m))-1)
#hsanames=mirnames[which(grepl("hsa",mirnames))]
#hsamature=maturenames[which(grepl("hsa",maturenames))]


prep_names=data.frame(hairpin_name=prep$Name)
prep_names$hairpin_seq=NA
prep_names$arm5p=prep$arm5p
prep_names$mature5p_seq=NA
prep_names$arm3p=prep$arm3p
prep_names$mature3p_seq=NA
options(warn=1)
for(i in 1:dim(prep_names)[1]){
  prep_names$hairpin_seq[i]=paste(s[grep(paste("^",prep_names$hairpin_name[i],"$",sep=""),mirnames)])
  prep_names$mature5p_seq[i]=paste(m[grep(paste("^",prep_names$arm5p[i],"$",sep=""),maturenames)])
  prep_names$mature3p_seq[i]=paste(m[grep(paste("^",prep_names$arm3p[i],"$",sep=""),maturenames)])
  
}

write.fasta(as.list(prep_names$hairpin_seq),prep_names$hairpin_name,"can.fa",open="w")
write.fasta(as.list(prep_names$mature5p_seq),prep_names$arm5p,"mat5p.fa",open="w")
write.fasta(as.list(prep_names$mature3p_seq),prep_names$arm3p,"mat3p.fa",open="w")
write.csv(prep_names,"prep_names.csv")
