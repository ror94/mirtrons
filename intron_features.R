intron_features <- function (input) {

patterns=matrix(,nrow=dim(input)[1],ncol=3)   
  
  for (i in 1:dim(input)[1]){
    patterns[i,1]=dim(str_locate_all(input$intron_seq[i],"TATA")[[1]])[1]/nchar(input$intron_seq[i])
    patterns[i,2]=dim(str_locate_all(input$intron_seq[i],"CAAT")[[1]])[1]/nchar(input$intron_seq[i])
    patterns[i,3]=dim(str_locate_all(input$intron_seq[i],"CG")[[1]])[1]/(str_count(input$intron_seq[i],"C")*str_count(input$intron_seq[i],"G")/
                                                                           nchar(input$intron_seq[i]))
    }


results=data.frame(intron_length=nchar(input$intron_seq))
results$intron5p_A=str_count(input$intron5p_seq,"A")*100/nchar(input$intron5p_seq)
results$intron5p_C=str_count(input$intron5p_seq,"C")*100/nchar(input$intron5p_seq)
results$intron5p_G=str_count(input$intron5p_seq,"G")*100/nchar(input$intron5p_seq)
results$intron5p_T=str_count(input$intron5p_seq,"T")*100/nchar(input$intron5p_seq)
results$intron3p_A=str_count(input$intron3p_seq,"A")*100/nchar(input$intron3p_seq)
results$intron3p_C=str_count(input$intron3p_seq,"C")*100/nchar(input$intron3p_seq)
results$intron3p_G=str_count(input$intron3p_seq,"G")*100/nchar(input$intron3p_seq)
results$intron3p_T=str_count(input$intron3p_seq,"T")*100/nchar(input$intron3p_seq)

results$intron5p_half5p_A=str_count(substr(input$intron5p_seq,1, ceiling(nchar(input$intron5p_seq)/2))
                                    ,"A")*100/nchar(substr(input$intron5p_seq,1, ceiling(nchar(input$intron5p_seq)/2)))
results$intron5p_half5p_C=str_count(substr(input$intron5p_seq,1, ceiling(nchar(input$intron5p_seq)/2))
                                    ,"C")*100/nchar(substr(input$intron5p_seq,1, ceiling(nchar(input$intron5p_seq)/2)))
results$intron5p_half5p_G=str_count(substr(input$intron5p_seq,1, ceiling(nchar(input$intron5p_seq)/2))
                                    ,"G")*100/nchar(substr(input$intron5p_seq,1, ceiling(nchar(input$intron5p_seq)/2)))
results$intron5p_half5p_T=str_count(substr(input$intron5p_seq,1, ceiling(nchar(input$intron5p_seq)/2))
                                    ,"T")*100/nchar(substr(input$intron5p_seq,1, ceiling(nchar(input$intron5p_seq)/2)))
results$intron3p_half5p_A=str_count(substr(input$intron3p_seq,1, ceiling(nchar(input$intron3p_seq)/2))
                                    ,"A")*100/nchar(substr(input$intron3p_seq,1, ceiling(nchar(input$intron3p_seq)/2)))
results$intron3p_half5p_C=str_count(substr(input$intron3p_seq,1, ceiling(nchar(input$intron3p_seq)/2))
                                    ,"C")*100/nchar(substr(input$intron3p_seq,1, ceiling(nchar(input$intron3p_seq)/2)))
results$intron3p_half5p_G=str_count(substr(input$intron3p_seq,1, ceiling(nchar(input$intron3p_seq)/2))
                                    ,"G")*100/nchar(substr(input$intron3p_seq,1, ceiling(nchar(input$intron3p_seq)/2)))
results$intron3p_half5p_T=str_count(substr(input$intron3p_seq,1, ceiling(nchar(input$intron3p_seq)/2))
                                    ,"T")*100/nchar(substr(input$intron3p_seq,1, ceiling(nchar(input$intron3p_seq)/2)))

results$intron5p_half3p_A=str_count(substr(input$intron5p_seq,ceiling(nchar(input$intron5p_seq)/2)+1,nchar(input$intron5p_seq)),
                                    "A")*100/nchar(substr(input$intron5p_seq,ceiling(nchar(input$intron5p_seq)/2)+1,nchar(input$intron5p_seq)))
results$intron5p_half3p_C=str_count(substr(input$intron5p_seq,ceiling(nchar(input$intron5p_seq)/2)+1,nchar(input$intron5p_seq)),
                                    "C")*100/nchar(substr(input$intron5p_seq,ceiling(nchar(input$intron5p_seq)/2)+1,nchar(input$intron5p_seq)))
results$intron5p_half3p_G=str_count(substr(input$intron5p_seq,ceiling(nchar(input$intron5p_seq)/2)+1,nchar(input$intron5p_seq)),
                                    "G")*100/nchar(substr(input$intron5p_seq,ceiling(nchar(input$intron5p_seq)/2)+1,nchar(input$intron5p_seq)))
results$intron5p_half3p_T=str_count(substr(input$intron5p_seq,ceiling(nchar(input$intron5p_seq)/2)+1,nchar(input$intron5p_seq)),
                                    "T")*100/nchar(input$intron5p_seq)
results$intron3p_half3p_A=str_count(substr(input$intron3p_seq,ceiling(nchar(input$intron3p_seq)/2)+1,nchar(input$intron3p_seq)),
                                    "A")*100/nchar(substr(input$intron3p_seq,ceiling(nchar(input$intron3p_seq)/2)+1,nchar(input$intron3p_seq)))
results$intron3p_half3p_C=str_count(substr(input$intron3p_seq,ceiling(nchar(input$intron3p_seq)/2)+1,nchar(input$intron3p_seq)),
                                    "C")*100/nchar(substr(input$intron3p_seq,ceiling(nchar(input$intron3p_seq)/2)+1,nchar(input$intron3p_seq)))
results$intron3p_half3p_G=str_count(substr(input$intron3p_seq,ceiling(nchar(input$intron3p_seq)/2)+1,nchar(input$intron3p_seq)),
                                    "G")*100/nchar(substr(input$intron3p_seq,ceiling(nchar(input$intron3p_seq)/2)+1,nchar(input$intron3p_seq)))
results$intron3p_half3p_T=str_count(substr(input$intron3p_seq,ceiling(nchar(input$intron3p_seq)/2)+1,nchar(input$intron3p_seq)),
                                    "T")*100/nchar(substr(input$intron3p_seq,ceiling(nchar(input$intron3p_seq)/2)+1,nchar(input$intron3p_seq)))
results$intron5p_fe=input$intron5p_fe
results$intron3p_fe=input$intron3p_fe
#results$introns_number=input$introns_number
#results$which_intron2=input$intron_number

results$tata=patterns[,1]
results$caat=patterns[,2]
results$CpG_ratio=patterns[,3]
#results$splice_donor=substr(input$intron_seq,1,3)
#results$splice_acc=substr(input$intron_seq,nchar(input$intron_seq)-2,nchar(input$intron_seq))
#results$tail5p=str_locate(input$intron_seq,results$Hairpin_seq)[,1]-1
#results$tail3p=nchar(input$intron_seq)-str_locate(input$intron_seq,results$Hairpin_seq)[,2]

return(results)
}
