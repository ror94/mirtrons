mirna_features <- function (input) {
  input$hairpin_seq=chartr("T","U",input$hairpin_seq)
  overhang_count=overhangcount(input)
  loops=loopscount(input)
  results=data_frame(hairpin_length = nchar(input$hairpin_seq))
  results = results %>% 
    tibble::rownames_to_column() %>% 
    mutate(mature5p_length=nchar(input$mature5p_seq),
    mature3p_length=nchar(input$mature3p_seq),
    mature5p_A=str_count(input$mature5p_seq,"A")*100/nchar(input$mature5p_seq),
    mature5p_C=str_count(input$mature5p_seq,"C")*100/nchar(input$mature5p_seq),
    mature5p_G=str_count(input$mature5p_seq,"G")*100/nchar(input$mature5p_seq),
    mature5p_U=str_count(input$mature5p_seq,"U")*100/nchar(input$mature5p_seq),
    mature3p_A=str_count(input$mature3p_seq,"A")*100/nchar(input$mature3p_seq),
    mature3p_C=str_count(input$mature3p_seq,"C")*100/nchar(input$mature3p_seq),
    mature3p_G=str_count(input$mature3p_seq,"G")*100/nchar(input$mature3p_seq),
    mature3p_U=str_count(input$mature3p_seq,"U")*100/nchar(input$mature3p_seq),
    interarm5p = str_locate(input$hairpin_seq, input$mature5p_seq)[,2]+1,
    interarm3p = str_locate(input$hairpin_seq, input$mature3p_seq)[,1]-1,
    interarm = substr(input$hairpin_seq,interarm5p,interarm3p),
    interarm_length = nchar(interarm),
    interarm_A=str_count(interarm,"A")*100/nchar(interarm),
    interarm_C=str_count(interarm,"C")*100/nchar(interarm),
    interarm_G=str_count(interarm,"G")*100/nchar(interarm),
    interarm_U=str_count(interarm,"U")*100/nchar(interarm),
    harpin_FE=input$fe/nchar(input$hairpin_seq),
    hairpin_A=str_count(input$hairpin_seq,"A")*100/nchar(input$hairpin_seq),
    hairpin_C=str_count(input$hairpin_seq,"C")*100/nchar(input$hairpin_seq),
    hairpin_G=str_count(input$hairpin_seq,"G")*100/nchar(input$hairpin_seq),
    hairpin_U=str_count(input$hairpin_seq,"U")*100/nchar(input$hairpin_seq),
    # mature5pposition[which(nchar(input$mature5p_seq)!=0)]=str_locate(input$hairpin_seq[which(nchar(input$mature5p_seq)!=0)],input$mature5p_seq[which(nchar(input$mature5p_seq)!=0)])[,1],
    # mature3pposition[which(nchar(input$mature3p_seq)!=0)]=hairpin_length[which(nchar(input$mature3p_seq)!=0)]-(str_locate(input$hairpin_seq[which(nchar(input$mature3p_seq)!=0)],input$mature3p_seq[which(nchar(input$mature3p_seq)!=0)])[,2])+1,
    overhang=overhang_count,
    small_loops=as.numeric(loops[,1]),
    large_loops=as.numeric(loops[,2]),
    t_loop_length = as.numeric(loops[,3]),
    class = input$mirna_class,
    hairpin_name = input$hairpin_name) %>%
    dplyr::select(-interarm)

return(results)
}
