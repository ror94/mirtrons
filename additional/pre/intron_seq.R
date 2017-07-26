input=read.csv("./Data/introns.csv", header=TRUE, stringsAsFactors = FALSE)
intron_seq=data.frame(matrix(,nrow=dim(input)[1],ncol=2))  
colnames(intron_seq)=c("intron5p_seq", "intron3p_seq")
for (i in 1:dim(input)[1]){
  x=nchar(input$intron_seq[i])
  seq=input$intron_seq[i]
  if (x<70){
    intron_seq$intron5p_seq[i]=seq
    intron_seq$intron3p_seq[i]=seq
  } else
    intron_seq$intron5p_seq[i]=substring(seq,1,70)
  intron_seq$intron3p_seq[i]=substring(seq,nchar(seq)-69,nchar(seq))
}

write.csv(intron_seq,"intronsiki.csv")