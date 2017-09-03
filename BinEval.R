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