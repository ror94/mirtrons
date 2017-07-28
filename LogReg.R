LogReg <-function(Z, folds){

set.seed(15)  
#X=BK_sample_crossVal(dim(Z)[1],itnumber)
LR=data.frame(Sensitivity=double(),Specificity=double(),F1=double(),AUC=double(), MCC=double())
RF=LR
LDA=LR
DT=LR
SVM=LR
NB=LR
#folds = generateCVRuns(Z$class, ntimes = 1, nfold = itnumber, stratified = TRUE)
for (i in 1:length(folds[[1]])){

  testing_data=Z[folds[[1]][[i]],]
  learning_data=Z[-folds[[1]][[i]],]
  check=testing_data$class
  testing_data$class=NULL
  
  
  
  # LOGISTIC REGRESSION
  model=glm(class~.,data=learning_data, family="binomial")
  model_pred_probs=predict(model, testing_data,type="response")
  model_results=rep(0,dim(testing_data)[1])
  model_results[model_pred_probs >(0.5)]=1
  df=data.frame(model_pred_probs, check)
  LR[i,]=BinEval(df)
  
  #RANDOM FOREST; default ntree=500
  rfm=randomForest(as.factor(class)~ ., learning_data, importance=TRUE)
  rfm_results=predict(rfm,testing_data)#, OOB = TRUE)
  df=data.frame(as.double(rfm_results)-1, check)
  table(rfm_results,check)
  RF[i,]=BinEval(df)
  
  #LDA
  ldam=lda(class~.,data=learning_data)
  ldam_results=predict(ldam,testing_data)$class
  df=data.frame(as.double(ldam_results)-1, check)
  LDA[i,]=BinEval(df)
  
  #Decision Tree
  tree=tree(as.factor(class)~.,learning_data)
  tree_results=predict(tree, testing_data, type="class")
  df=data.frame(as.double(tree_results)-1, check)
  #pruned_tree=prune.misclass(tree, best=8) #pruning
  #pruned_tree_results=predict(pruned_tree, testing_data, type="class")
  #df=data.frame(as.double(pruned_tree_results)-1, check)
  DT[i,]=BinEval(df)
  #print(DT[i,])
  
  #SVM
  svm_model=svm(as.factor(class)~.,data=learning_data)
  svm_results=predict(svm_model,testing_data,type="class")
  df=data.frame(as.double(svm_results)-1, check)
  #tuned=tune(svm, as.factor(class)~.,data=learning_data, ranges = list(epsilon = seq(0,0.2,0.01), cost = c(0.001,0.01,.1,1,10,100)))
  SVM[i,]=BinEval(df)
  
  #Naive Bayes
  nb=naiveBayes(as.factor(class)~.,data=learning_data)
  nb_results=predict(nb, testing_data,type="class")
  df=data.frame(as.double(nb_results)-1, check)
  NB[i,]=BinEval(df)
  
  
}
LR=colSums(LR)/itnumber
RF=colSums(RF)/itnumber
LDA=colSums(LDA)/itnumber
DT=colSums(DT)/itnumber
SVM=colSums(SVM)/itnumber
NB=colSums(NB)/itnumber
Results=rbind(LR,RF,LDA,DT,SVM,NB)
rownames(Results)=c("Logistic Regression", "Random Forest", "Linear Discriminant Analysis","Decision Tree", "Support Vector machines","Naive Bayes")
results=list(Results,rfm,svm_model,folds)
return(results)
}
