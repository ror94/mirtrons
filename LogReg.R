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
