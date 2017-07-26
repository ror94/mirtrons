setwd("/home/ror/Desktop/Praca Inżynierska/NOWE")
mirna_coords=read.delim("hsa.gff3", header=F, comment.char="#")
mirna_names=read.table("mirnalist.txt")

coordinates=data.frame(matrix(,nrow=dim(mir_names)[1], ncol=5))
colnames(coordinates)=c("miRNA","Chromosome","Start","End","Strand")
  for (i in 1:dim(mirna_names)[1]){
    k=which(!is.na(str_locate(mirna_coords[,9], toString(mirna_names[i,1]))))
    k=k[1]
    coordinates$miRNA[i]=toString(mirna_names[i,1])
    coordinates$Chromosome[i]=toString(mirna_coords[k,1])
    coordinates$Start[i]=toString(mirna_coords[k,4])
    coordinates$End[i]=toString(mirna_coords[k,5])
    coordinates$Strand[i]=toString(mirna_coords[k,7])
    
  }

write.csv(coordinates,"coordinates.csv")
