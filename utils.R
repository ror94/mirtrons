

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
loopscount <- function(input){
  
  all_loops=matrix(,dim(input)[1],3)
  input$hairpin_seq=chartr("T","U",input$hairpin_seq)
  errors=c()
  
  for (k in 1:dim(input)[1]){
    tryCatch({
      db=input$dotbracket[k]
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
  colnames(all_loops) = c("small_loops", "large_loops", "terminal_loop_length")
  return(all_loops)
}
LogReg <-function(Z, folds, models = "all"){
  
  #X=BK_sample_crossVal(dim(Z)[1],itnumber)
  LR=data_frame(Sensitivity=double(),Specificity=double(),F1=double(),AUC=double(), MCC=double())
  RF=LR
  LDA=LR
  DT=LR
  SVM=LR
  NB=LR
  rfm = c()
  #folds = generateCVRuns(Z$class, ntimes = 1, nfold = itnumber, stratified = TRUE)
  for (i in 1:length(folds[[1]])){
    testing_data=Z[folds[[1]][[i]],]
    learning_data=Z[-folds[[1]][[i]],]
    check=testing_data$class
    testing_data$class=NULL
    
    if (models == "LR" | models == "all"){
      # LOGISTIC REGRESSION
      model=glm(class~.,data=learning_data, family="binomial")
      model_pred_probs=predict(model, testing_data,type="response")
      df=data.frame(model_pred_probs, check)
      LR[i,]=BinEval(df)
    }
    if (models == "RF" | models == "all"){
      #RANDOM FOREST; default ntree=500
      rfm=randomForest(as.factor(class)~ ., learning_data, importance=TRUE)
      rfm_results=predict(rfm,testing_data)#, OOB = TRUE)
      df=data.frame(as.double(rfm_results)-1, check)
      table(rfm_results,check)
      RF[i,]=BinEval(df)
    }
    if (models == "LDA" | models == "all"){
      #LDA
      ldam=lda(class~.,data=learning_data)
      ldam_results=predict(ldam,testing_data)$class
      df=data.frame(as.double(ldam_results)-1, check)
      LDA[i,]=BinEval(df)
    }
    if (models == "DT" | models == "all"){
      #Decision Tree
      tree=tree(as.factor(class)~.,learning_data)
      tree_results=predict(tree, testing_data, type="class")
      df=data.frame(as.double(tree_results)-1, check)
      #pruned_tree=prune.misclass(tree, best=8) #pruning
      #pruned_tree_results=predict(pruned_tree, testing_data, type="class")
      #df=data.frame(as.double(pruned_tree_results)-1, check)
      DT[i,]=BinEval(df)
      #print(DT[i,])
    }
    if (models == "SVM" | models == "all"){
      #SVM
      svm_model=svm(as.factor(class)~.,data=learning_data)
      svm_results=predict(svm_model,testing_data,type="class")
      df=data.frame(as.double(svm_results)-1, check)
      #tuned=tune(svm, as.factor(class)~.,data=learning_data, ranges = list(epsilon = seq(0,0.2,0.01), cost = c(0.001,0.01,.1,1,10,100)))
      SVM[i,]=BinEval(df)
    }
    if (models == "NB" | models =="all"){
      #Naive Bayes
      nb=naiveBayes(as.factor(class)~.,data=learning_data)
      nb_results=predict(nb, testing_data,type="class")
      df=data.frame(as.double(nb_results)-1, check)
      NB[i,]=BinEval(df)
    }
    
  }
  LR=colSums(LR)/itnumber
  RF=colSums(RF)/itnumber
  LDA=colSums(LDA)/itnumber
  DT=colSums(DT)/itnumber
  SVM=colSums(SVM)/itnumber
  NB=colSums(NB)/itnumber
  Results = bind_rows(LR,RF,LDA,DT,SVM,NB)%>%
    mutate(Method = c("Logistic Regression", "Random Forest", 
                      "Linear Discriminant Analysis","Decision Tree", "Support Vector Machines","Naive Bayes"))   %>%
    arrange(desc(MCC)) %>%
    .[c(6,1,2,4,3,5)]
  results=list(results = Results,rf = rfm,svm = svm_model,folds = folds)
  return(results)
}
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
BinEval <-function (df){
  pred=ROCR::prediction(df[1],df[2])
  auc=ROCR::performance(pred, "auc")
  sens = ROCR::performance(pred, "sens")
  spec = ROCR::performance(pred, "spec")
  roc = ROCR::performance(pred, "tpr", "fpr")
  #lift = ROCR::performance(pred, "lift", "tpr")
  
  
  model_results=rep("Canonical",dim(df)[1])
  model_results[which(df[,1] >(0.5))]="Mirtron"
  x=table(model_results,df[,2])
  TN = sum(model_results == "Canonical" & df[,2] == "Canonical")
  TP = sum(model_results == "Mirtron" & df[,2] == "Mirtron")
  FP = sum(model_results == "Mirtron" & df[,2] == "Canonical")
  FN = sum(model_results == "Canonical" & df[,2] == "Mirtron")
  #mean(model_results != check)
  Sens=TP/(TP+FN)
  Spec=TN/(TN+FP)
  Acc=(TP+TN)/sum(x)
  F1=2*TP/(2*TP+FP+FN)
  MCC=(TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
  
  x = confusionMatrix(model_results, df[[2]], positive = "Mirtron")
  Sens = x$byClass["Sensitivity"]
  Spec = x$byClass["Specificity"]
  Acc = x$byClass["Precision"]
  F1 = x$byClass["F1"]
  
  Results=data.frame(Sensitivity=Sens,Specificity=Spec,F1=F1,AUC=auc@y.values[[1]],Mathew=MCC)
  
  return(Results)
}
grbiplot <- function(PC, class, ...){
  plot_data = data.frame(PC1 =PC$x[,1], PC2 = PC$x[,2], class = class)
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(plot_data$PC2) - min(plot_data$PC2)/(max(datapc$PC2)-min(datapc$PC2))),
    (max(plot_data$PC1) - min(plot_data$PC1)/(max(datapc$PC1)-min(datapc$PC1)))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get("PC1")),
                      v2 = .7 * mult * (get("PC2"))
  )
  datapc$angle <- with(datapc, (180/pi) * atan(v2/v1))
  datapc$hjust = with(datapc, (1 - 1 * sign(v1))/2)
  
  g = ggplot(plot_data)+
    geom_point(aes(x = PC1, y = PC2, color = class), size = 2)+
    geom_text_repel(data=datapc, aes(x=v1, y=v2, label=varnames),#, angle = angle),
                    color=muted("red"), size = 5, alpha = 1, segment.size = 0.3)+
    geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), 
                 arrow=arrow(length=unit(0.2,"cm")), alpha=1, color=muted("red"), size = 0.6) +
    theme(legend.justification=c(1,0), legend.position=c(1,0))+ ...
  print(g)
}
stepwise <- function(ml_data, ...){
  preds = ml_data %>% select(-class) %>% names
  stepwise = c()
  f1s = c()
  models_table = list()
  SVM = data_frame(Sensitivity=double(),Specificity=double(),F1=double(),AUC=double(), MCC=double())
  #models = "SVM"
  #stepwise_svm = function(folds, models){
  for (i in 1:length(preds)){
    f1 = 0
    for (j in setdiff(preds, stepwise)){
      Z = ml_data %>% select(c(stepwise, j), class) 
      xs=LogReg(Z, folds, models = "SVM")
      F1 = xs[[1]] %>% filter(Method == "Support Vector Machines") %>% select(F1)
      if (!is.na(F1) & F1 > f1 ) {
        best_model <- xs[[3]]
        f1 <- as.numeric(F1)
        best_variable <- j
      }
    }
    stepwise = append(stepwise, best_variable)
    f1s = append(f1s, f1)
    models_table[[i]] = best_model
    cat(paste0(i,". ",best_variable, ", F1 = ", f1, "\n"))
  }
  stepwise_results = data_frame(feature = stepwise, "F1" = f1s) %>% rownames_to_column()
  return(stepwise_results)
  
  
  
}
