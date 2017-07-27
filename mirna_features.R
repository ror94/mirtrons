mirna_features <- function (input,random) {
  input$hairpin_seq=chartr("T","U",input$hairpin_seq)
  results=data.frame(hairpin_length=nchar(input$hairpin_seq))
  results$mature5p_length=nchar(input$mature5p_seq)
  results$mature3p_length=nchar(input$mature3p_seq)
  results$mature5p_A=str_count(input$mature5p_seq,"A")*100/nchar(input$mature5p_seq)
  results$mature5p_C=str_count(input$mature5p_seq,"C")*100/nchar(input$mature5p_seq)
  results$mature5p_G=str_count(input$mature5p_seq,"G")*100/nchar(input$mature5p_seq)
  results$mature5p_U=str_count(input$mature5p_seq,"U")*100/nchar(input$mature5p_seq)
  results$mature3p_A=str_count(input$mature3p_seq,"A")*100/nchar(input$mature3p_seq)
  results$mature3p_C=str_count(input$mature3p_seq,"C")*100/nchar(input$mature3p_seq)
  results$mature3p_G=str_count(input$mature3p_seq,"G")*100/nchar(input$mature3p_seq)
  results$mature3p_U=str_count(input$mature3p_seq,"U")*100/nchar(input$mature3p_seq)
  results$mature3p_A[is.nan(results$mature3p_A)]=0
  results$mature3p_C[is.nan(results$mature3p_C)]=0
  results$mature3p_G[is.nan(results$mature3p_G)]=0
  results$mature3p_U[is.nan(results$mature3p_U)]=0
  results$mature5p_A[is.nan(results$mature5p_A)]=0
  results$mature5p_C[is.nan(results$mature5p_C)]=0
  results$mature5p_G[is.nan(results$mature5p_G)]=0
  results$mature5p_U[is.nan(results$mature5p_U)]=0
  interarm5p = str_locate(input$hairpin_seq, input$mature5p_seq)[,2]+1
  interarm3p = str_locate(input$hairpin_seq, input$mature3p_seq)[,1]-1
  interarm = substr(input$hairpin_seq,interarm5p,interarm3p)
  results$interarm_length = nchar(interarm)
  results$interarm_A=str_count(interarm,"A")*100/nchar(interarm)
  results$interarm_C=str_count(interarm,"C")*100/nchar(interarm)
  results$interarm_G=str_count(interarm,"G")*100/nchar(interarm)
  results$interarm_U=str_count(interarm,"U")*100/nchar(interarm)
  results$interarm_A[is.nan(results$interarm_A)]=0
  results$interarm_C[is.nan(results$interarm_C)]=0
  results$interarm_G[is.nan(results$interarm_G)]=0
  results$interarm_U[is.nan(results$interarm_U)]=0
  if (random==FALSE){
    overhang=overhangcount(input)
    loops=loopscount(input)
    results$Harpin_FE=as.numeric(input$fe/nchar(input$hairpin_seq))
    #results$Harpin_FE=as.numeric(input$fe)
    results$hairpin_A=str_count(input$hairpin_seq,"A")*100/nchar(input$hairpin_seq)
    results$hairpin_C=str_count(input$hairpin_seq,"C")*100/nchar(input$hairpin_seq)
    results$hairpin_G=str_count(input$hairpin_seq,"G")*100/nchar(input$hairpin_seq)
    results$hairpin_U=str_count(input$hairpin_seq,"U")*100/nchar(input$hairpin_seq)
    results$mature5pposition[which(nchar(input$mature5p_seq)!=0)]=str_locate(input$hairpin_seq[which(nchar(input$mature5p_seq)!=0)],input$mature5p_seq[which(nchar(input$mature5p_seq)!=0)])[,1]
    results$mature3pposition[which(nchar(input$mature3p_seq)!=0)]=results$hairpin_length[which(nchar(input$mature3p_seq)!=0)]-(str_locate(input$hairpin_seq[which(nchar(input$mature3p_seq)!=0)],input$mature3p_seq[which(nchar(input$mature3p_seq)!=0)])[,2])+1
    results$overhangnumber=overhang
    results$small_loops=as.numeric(loops[,1])
    results$large_loops=as.numeric(loops[,2])
    results$t_loop_length = as.numeric(loops[,3])
  }
  else{}
  
return(results)
}