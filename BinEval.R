BinEval <-function (df){
  pred=ROCR::prediction(df[1],df[2])
  auc=ROCR::performance(pred, "auc")
  
  
  model_results=rep(0,dim(df)[1])
  model_results[which(df[,1] >(0.5))]=1
  x=table(model_results,df[,2])
  TN=x[1]
  TP=x[4]
  FP=x[2]
  FN=x[3]
  #mean(model_results != check)
  Sens=TP/(TP+FN)
  Spec=TN/(TN+FP)
  Acc=(TP+TN)/sum(x)
  F1=2*TP/(2*TP+FP+FN)
  MCC=(TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
  Results=data.frame(Sensitivity=Sens,Specificity=Spec,F1=F1,AUC=auc@y.values[[1]],Mathew=MCC)
  
  return(Results)
}