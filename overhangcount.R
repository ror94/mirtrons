overhangcount <- function(input){
  
  overhang=matrix(,dim(input)[1],1)
  input$hairpin_seq=chartr("T","U",input$hairpin_seq)
  
  for (i in 1:dim(input)[1]){
    
    mir_seq=input$hairpin_seq[i]
    db=input$dotbracket[i]
    mir5p=input$mature5p_seq[i]
    mir3p=input$mature3p_seq[i]
    
    if (nchar(mir3p)==0){
      overhang[i]=""
     } else if (nchar(mir5p)==0){
      overhang[i]=""
    } else { 
      mir5pstart=str_locate(mir_seq,mir5p)[1]
      mir5pend=str_locate(mir_seq,mir5p)[2]
      mir3pstart=str_locate(mir_seq,mir3p)[1]
      mir3pend=str_locate(mir_seq,mir3p)[2]
      x=makeCt(db,mir_seq)
      y=matrix(,2,dim(x)[1])
      y[1,]=x$pos2
      y[2,]=x$bound
      mir5pct=y[,mir5pstart:mir5pend]
      mir3pct=y[,mir3pstart:mir3pend]
      mir3pct=mir3pct[2:1,ncol(mir3pct):1]
      overhang[i]=which(mir5pct[1,]==intersect(mir5pct[1,],mir3pct[1,])[1])-which(mir3pct[1,]==intersect(mir5pct[1,],mir3pct[1,])[1]) #find first common elemen
     #+ is overhang 3p
     #- is ovehang 5p
    }
    
    
  }  
  
  return(as.numeric(overhang))
}

