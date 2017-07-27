loopscount <- function(input){

all_loops=matrix(,dim(input)[1],3)
input$hairpin_seq=chartr("T","U",input$hairpin_seq)
errors=c()
  
    for (k in 1:dim(input)[1]){
      tryCatch({
        db=input$db[k]
        seq=input$hairpin_seq[k]
        #db="(((.((((..(((....))).)))..))))"
        #seq="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        
        ct=makeCt(db,seq)
        
        
        last5p=tail(which((ct$pos<ct$bound & ct$bound!=0)==T),n=1)
        first3p=which((ct$pos>ct$bound & ct$bound!=0)==T & ct$pos > last5p)[1]
        
        x=str_locate_all(db,"\\(")[[1]][,1]
        y=str_locate_all(db,"\\)")[[1]][,1]
        if (min(y)<(0.4*nchar(db))){
        x=x[which(y==ct$bound[max(x)]):length(x)] #cuts when loops at the beginning
        y=y[which(y==ct$bound[max(x)]):length(y)] #cuts when loops at the beginning
        } else if(max(x)>(0.6*nchar(db))) {
        x=x[1:which(x==which(ct$bound==min(y)))] #cuts when loops at the beginning
        y=y[1:which(x==which(ct$bound==min(y)))] #cuts when loops at the beginning 
          
        }
        arm5p=with(ct,rbind(pos[min(x):max(x)],bound[min(x):max(x)]))
        arm3p=with(ct,rbind(rev(pos[min(y):max(y)]),rev(bound[min(y):max(y)])))
        arm3p=apply(arm3p,2,rev)
        diffs=setdiff(arm5p[1,],arm3p[1,])
        if (length(diffs)>0){
            for (i in 1:length(diffs)){
              pre=arm3p[,1:which(arm3p[1,]==(diffs[i]-1))[1]]
              post=arm3p[,(which(arm3p[1,]==(diffs[i])-1)+1)[1]:length(arm3p[1,])]
              cp=matrix(arm5p[,which(arm5p[1,]==diffs[i])],nrow=2)
              arm3p=cbind(cbind(pre,cp),post)
            }
            hairpin=arm3p
            log_hairpin=hairpin[1,] & hairpin[2,]
            log_hairpin=as.numeric(!log_hairpin)
            log_hairpin=paste(log_hairpin,sep="",collapse="")
            loops_ind=str_locate_all(log_hairpin,"1+")[[1]]
            loops_number=length(loops_ind[,1])
            loop_size=c()
            for (i in 1:loops_number){
              piece=hairpin[,loops_ind[i,1]:loops_ind[i,2]]
             log_piece=matrix(hairpin[,loops_ind[i,1]:loops_ind[i,2]]==0,nrow = 2)
             loop_size[i]=max(apply(log_piece,1,sum))
            }
            small_loops=sum(loop_size<4)
            large_loops=sum(loop_size>=4)
            all_loops[k,]=c(small_loops,large_loops, "none")
        } else {
            all_loops[k,]=c(0,0, "none")
            }
        all_loops[k,3] = first3p - last5p - 1
  }, error = function(err){
    print(k)
  }) 
    }
return(all_loops)
}
