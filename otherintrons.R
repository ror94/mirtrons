introns=readDNAStringSet("introns.fa")
#x=which(nchar(paste(introns))<200)
x=1:length(introns)
all_introns=data.frame(matrix(nrow=length(x),ncol=5))
colnames(all_introns)=c("tx_name", "intron_seq", "intron_name", "intron_number","introns_number")
for (k in 1:length(x)){
  i=x[k]
  all_introns$tx_name[k]=substr(introns@ranges@NAMES[i],str_locate(introns@ranges@NAMES[i],"Gene_")[2]+1,str_locate(introns@ranges@NAMES[i],"Gene_")[2]+8)
  all_introns$intron_seq[k]=paste(introns)[i]
  all_introns$introns_number[k]=length(which(!is.na(str_locate(introns@ranges@NAMES,all_introns$tx_name[k]))[,1]))
  all_introns$intron_name[k]=introns@ranges@NAMES[i]
  strand=substr(introns@ranges@NAMES[k],str_locate(introns@ranges@NAMES[i],"strand=")[2]+1,str_locate(introns@ranges@NAMES[i],"strand=")[2]+1)
  if (strand=="+") {
    all_introns$intron_number[k]=as.numeric(substr(all_introns$intron_name[k],27,28))+1
  }
  else {
    all_introns$intron_number[k]=all_introns$introns_number[k]-as.numeric(substr(all_introns$intron_name[k],27,28))
  }
}

all_other_introns=all_introns[ !(all_introns$intron_name %in% introns_df$intron_name), ]
write.csv(all_other_introns, "done.csv")
