classes = ml_data$class
mcc <- function(x, lev = NULL, model = NULL) {
  #f1_val <- F1_Score(y_pred = data$pred, y_true = data$obs, positive = lev[1])
  conf_mat = ConfusionMatrix(y_pred = x$pred, y_true = x$obs)
  TP = as.numeric(conf_mat[4])
  TN = as.numeric(conf_mat[1])
  FP = as.numeric(conf_mat[3])
  FN = as.numeric(conf_mat[2])
  c(F1 = 2*TP/(2*TP+FP+FN))
  #c(TPR = TP/(TP+FN))
  #(TNR = TN / (TN+FP))
  c(MCC = (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
}
caretFuncs$summary <- mcc
rfe(pca_data,classes,
                  sizes = c(1:21),
                  rfeControl = rfeControl(functions = caretFuncs, method = "cv",
                                          number = 5),
                  ## pass options to train()
                  method = "svmRadial", metrics = "F1")



plot(svmProfile)
predictors(svmProfile)
plot(svmProfile, type=c("g", "o"))
