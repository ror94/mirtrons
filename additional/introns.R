introns <- function(data){
introns=readDNAStringSet("./Fasta/introns.fa")
introns_df=data.frame(matrix(nrow=dim(data)[1],ncol=5))
colnames(introns_df)=c("tx_name", "intron_seq", "intron_name", "intron_number","introns_number")

  for (i in 1:nrow(data)){
    index=which(!is.na(str_locate(paste(introns),data$hairpin_seq[i])))[1] #finds given hairpin sequence
    introns_df$tx_name[i]=data$ucsc_id[i]
    introns_df$intron_seq[i]=paste(introns)[index]
    introns_df$introns_number[i]=length(which(!is.na(str_locate(introns@ranges@NAMES,data$ucsc_id[i]))[,1]))
    introns_df$intron_name[i]=introns@ranges@NAMES[i]
      if (data$strand[i]=="+") {
        introns_df$intron_number[i]=as.numeric(substr(introns_df$intron_name[i],27,28))+1
      }
      else {
        introns_df$intron_number[i]=introns_df$introns_number[i]-as.numeric(substr(introns_df$intron_name[i],27,28))
      }
}

return(introns_df)
}