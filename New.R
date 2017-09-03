rm(list = ls())
if (Sys.info()[[1]]=="Windows"){
  setwd("D:/BACKUP/Workspace/R/Mirtrons")
} else{
  setwd("~/Workspace/R/Mirtrons")
}

# PACKAGES ####
library(magrittr)
library(stringr) 
library("RRNA") # makect
library(agricolae) # tukey
library(ggbiplot)
library(devtools)
library(TunePareto) #CV runs
library(randomForest)
library(MASS) #lda
library(tree) #tree
library(e1071) #svm
library(ROCR) #evaluation
library(tibble)
library(dplyr)
library(tidyr)
library(Boruta)
library(caret) #RFE
library(MLmetrics) #cnfusion matrxi
library(lattice)
library(latticeExtra)
install_github("vqv/ggbiplot")

source("mirna_features.R")
source("overhangcount.R")
source("loopscount.R")
source("LogReg.R")
source("BinEval.R")

## DATA PREPERATION ####
canonical_data = read.csv("./Data/data.csv", header=TRUE, stringsAsFactors = FALSE)
canonical_data = tbl_df(canonical_data)
canonical_data %>% 
  tibble::rownames_to_column() %>% 
  mutate(mirna_class = factor(class, labels = c('Canonical', 'Mirtron'))) %>%
  dplyr::select(-class) %>%
  filter(!(as.numeric(rowname) %in% c(380,702,720,813,889,928))) %>%
  mirna_features -> mirna_data

test_mirna_data=read.csv("./Data/testdata.csv", header=TRUE, stringsAsFactors = FALSE)
test_mirna_data %>% 
  mutate(mirna_class = factor(class, labels = 'Mirtron')) %>%
  tibble::rownames_to_column() %>% 
  filter(!(as.numeric(rowname) %in% c(1,22,103,139,151,164,165,182,202))) %>% 
  mirna_features -> test_data

## PLOTS AND STATISTICS ####
source('mirnaplots.R')
#mirtronplots_mirna=mirnaplots(mirtron_mirna)
#canonicalplots_mirna=mirnaplots(canonical_mirna)
#testplots_mirna=mirnaplots(test_mirna)

mirtron_mature5p_G = mirna_data %>% filter(class == "Mirtron") %>% select(mature5p_G)
canonical_mature5p_G = mirna_data %>% filter(class == "Canonical") %>% select(mature5p_G)
mirtron_mature3p_C = mirna_data %>% filter(class == "Mirtron") %>% select(mature3p_C)

plot(ecdf(pull(mirtron_mature5p_G)))
plot(ecdf(pull(canonical_mature5p_G)))

plot(ecdf(pull(mirtron_mature5p_G)), verticals=TRUE, do.points=FALSE)
plot(ecdf(pull(canonical_mature5p_G)), verticals=TRUE, do.points=FALSE, add=TRUE, col='brown')

P.values = data_frame()
sign=c('Not Significant','Significant')

tukey = function(x, classes) {
  sign=c('Not Significant','Significant')
  tukey_df=data_frame(input = x, classes = classes)
  means=tapply(tukey_df$input,tukey_df$classes, mean)
  ajuste <- lm( tukey_df$input ~ tukey_df$classes)
  tukey.alpha=0.005
  h=HSD.test(ajuste, 'tukey_df$classes',alpha=tukey.alpha)
  tukey = sign[length(levels(h$groups$M))]
}
classes = as.numeric(mirna_data$class == "Mirtron")
mirna_data %>% 
  summarise_if(is.numeric, funs("T - test" = t.test(.[class == "Mirtron"], .[class == "Canonical"])$p.value,
                                "KS - test" = ks.test(.[class == "Mirtron"], .[class == "Canonical"])$p.value,
                                "Tukey - test" = tukey(., classes))) %>%
  gather %>%
  transmute(name = sub("\\_[^\\_]*$", "", key), test = sub(".*_", "", key), value = value) %>%
  spread(test, value) -> test_results

cat("\nStatistical tests\n")
print(test_results)

## Principal Component Analysis ####

ml_data = mirna_data %>% 
  dplyr::select(-c(contains("position"), contains("_U"), rowname, hairpin_name, interarm3p, interarm5p))
pca_data = dplyr::select(ml_data, -class)
pca=prcomp(pca_data, retx=TRUE, center=TRUE, scale=TRUE)
labels=factor(mirna_data$class)
g = ggbiplot(pca, obs.scale = 1, var.scale = 1, 
              groups = labels, ellipse = F, 
              circle = F) + scale_color_discrete(name = '') + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g) #mirtron vs canonical


ml_test_data = test_data %>%
  dplyr::select(-c(contains("position"), contains("_U"), rowname, hairpin_name, interarm3p, interarm5p))
pca_test_data = dplyr::select(ml_test_data, -class)
pca_test_data = bind_rows(pca_data, pca_test_data)
pca2 = pca
pca2$x = scale(pca_test_data, pca$center, pca$scale) %*% pca$rotation
labels = factor(c(as.vector(mirna_data$class), rep('test',nrow(pca_test_data) - nrow(pca_data))))
g2 <- ggbiplot(pca2, obs.scale = 1, var.scale = 1, 
              groups = labels, ellipse = F, 
              circle = F) + scale_color_discrete(name = '') + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g2) #added test

# 3d PCA
library(rgl)
colors = replace(as.vector(ml_data$class), with(ml_data, which(class == "Mirtron")), "blue")
colors = replace(colors,which(colors != "blue"), "red")
coords <- NULL
for (i in 1:nrow(pca$rotation)) {
  coords <- rbind(coords, rbind(c(0,0,0),pca$rotation[i,1:3]))
}
pca3 = plot3d(pca$x[,1:3], col=colors, size = 5)
text3d(pca$rotation[,1:3]*10, texts=rownames(pca$rotation), col="red")
lines3d(coords*10, col="red", lwd=2)

## CLASSIFICATION ####
source("LogReg.R")
source("BinEval.R")
itnumber = 5
set.seed(27)
folds = generateCVRuns(ml_data$class, ntimes = 1, nfold = itnumber, stratified = TRUE)
cat("\n1x5-fold Cross-validation classification\n")
x=LogReg(ml_data, folds)
print(x$results)
pred=predict(x$svm,ml_test_data) #predict on svm model
print(table(pred))
print(mean(as.double(pred)-1))

## SINGLE FEATURE ####

Singlef=data_frame()
for (i in 1:(dim(ml_data)[2]-1)){
  singlef=LogReg(ml_data[,c(i,ncol(ml_data))],folds)
  Singlef=bind_rows(Singlef,singlef$results %>% filter(grepl('Support', Method)))
}
Singlef = Singlef %>%
  mutate(Method = colnames(pca_data)) %>%
  arrange(desc(MCC))
print(Singlef)

## BORUTA ####
Bor = Boruta(ml_data[,-(ncol(ml_data))],ml_data$class, getImp = getImpRfZ)
plot(Bor)
labs = rev(names(colMeans(Bor$ImpHistory)))
#text(cex=1, x=x-.25, y=-1.25, labs, xpd=TRUE, srt=45)
Boruta_results = data_frame(Feature = colnames(Bor$ImpHistory), Means = colMeans(Bor$ImpHistory)) %>%
  arrange(desc(Means))

## STEPWISE SVM ####
preds = ml_data %>% select(-class) %>% names
stepwise = c()
f1s = c()
models_table = list()
SVM = data_frame(Sensitivity=double(),Specificity=double(),F1=double(),AUC=double(), MCC=double())
for (i in 1:length(preds)){
  f1 = 0
  for (j in setdiff(preds, stepwise)){
    Z = ml_data %>% select(c(stepwise, j), class) 
    x=LogReg(Z, folds, models = "SVM")
    F1 = x[[1]] %>% filter(Method == "Support Vector Machines") %>% select(F1)
    if (F1 > f1) {
      best_model <- x[[3]]
      f1 <- as.numeric(F1)
      best_variable <- j
    }
  }
  stepwise = append(stepwise, best_variable)
  f1s = append(f1s, f1)
  models_table[[i]] = best_model
  cat(i,". ",best_variable, ", F1 = ", f1, "\n")
}
stepwise_results = data_frame(feature = stepwise, "F1" = f1s) %>% rownames_to_column()


g = ggplot(stepwise_results, aes(x = feature, F1)) + 
  geom_point() +
  scale_x_discrete(limits = stepwise_results$feature) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(g)

## VERIFICATION ####
no_best_features = which(stepwise_results$F1 == max(stepwise_results$F1))
final_data = ml_data %>% select(stepwise_results$feature[1:no_best_features], class)
x=LogReg(final_data, folds)
print(x$results)

## RFE ####
classes = ml_data$class
library(caret) # RFE
mcc <- function(x, lev = NULL, model = NULL) {
  #f1_val <- F1_Score(y_pred = data$pred, y_true = data$obs, positive = lev[1])
  conf_mat = ConfusionMatrix(y_pred = x$pred, y_true = x$obs)
  TP = as.numeric(conf_mat[4])
  TN = as.numeric(conf_mat[1])
  FP = as.numeric(conf_mat[3])
  FN = as.numeric(conf_mat[2])
  c(F1 = 2*TP/(2*TP+FP+FN))
  #c(TPR = TP/(TP+FN))
  #c(TNR = TN / (TN+FP))
  #c(MCC = (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
}
caretFuncs$summary <- mcc

ctrl <- rfeControl(functions = caretFuncs,
                   method = "cv",
                   verbose = FALSE,
                   index = folds$`Run  1`)
svmProfile <-rfe(as.data.frame(pca_data), classes, sizes=c(1:21), rfeControl = ctrl, method = "svmRadial", metrics = "F1")

plot(svmProfile)
predictors(svmProfile)
plot(svmProfile, type=c("g", "o"))

## Print results ####
print(test_results, n = Inf)
print(x$results, n = Inf)
print(Singlef, n = Inf)
print(stepwise_results, n = Inf)
print(Boruta_results, n = Inf)
print(bind_cols(Boruta_results %>% filter(!grepl("shadow",Feature)), stepwise_results), n = Inf)

