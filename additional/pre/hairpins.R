
  mirna_names=read.table("canonicalmirna")
  hairpins= readRNAStringSet("hairpin.fa")
  hairpins_table=data.frame(matrix(,nrow=dim(mirna_names)[1], ncol=2))
  colnames(hairpins_table)=c("hairpin_name","hairpin_seq")
  
  for (i in 1:dim(mirna_names)[1]){
    name=toString(mirna_names[i,1])
    name=paste(name, "")
    index=which(!is.na(str_locate(names(hairpins),name)[,1]), TRUE)
    hairpins_table[i,1]=name
    hairpins_table[i,2]=paste(hairpins)[index]
  }
  
  hairpins_table$mir_name=gsub("-1 ","",hairpins_table$hairpin_name)
  hairpins_table$mir_name=gsub("-2 ","",hairpins_table$mir_name)
  hairpins_table$mir_name=gsub("-3 ","",hairpins_table$mir_name)
  hairpins_table$mir_name=gsub(" ","",hairpins_table$mir_name)
  hairpins_table$mir_name=chartr("","-",hairpins_table$mir_name)
  hairpins_table$mir_name=chartr("r","R",hairpins_table$mir_name)
  
write.csv(hairpins_table, "hairpindata.csv", row.names=F)

