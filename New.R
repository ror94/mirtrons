rm(list = ls())
if (Sys.info()[[1]]=="Windows"){
  setwd("D:/BACKUP/Workspace/R/Mirtrons")
} else{
  setwd("~/Workspace/R/Mirtrons")
}

## PACKAGES ####
#source("https://bioconductor.org/biocLite.R")
#biocLite("RRNA")
library(devtools)
install_github("vqv/ggbiplot")
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
library(tidyverse)
library(Boruta)
library(caret) #RFE
library(lattice)
library(grid)
library(gridExtra)
#library(latticeExtra)
source("mirna_features.R")
source("overhangcount.R")
source("loopscount.R")
source("LogReg.R")
source("BinEval.R")
source("multiplot.R")


## DATA PREPERATION ####
canonical_data = read.csv("./Data/data.csv", header=TRUE, stringsAsFactors = FALSE)
canonical_data = tbl_df(canonical_data)
canonical_data %>% 
  tibble::rownames_to_column() %>% 
  mutate(mirna_class = factor(class, labels = c('Canonical', 'Mirtron'))) %>%
  dplyr::select(-class) %>%
  filter(!(as.numeric(rowname) %in% c(380,702,720,813,889,928))) %>%
  mirna_features %>%
  select(-c(interarm3p, interarm5p)) -> mirna_data

test_mirna_data=read.csv("./Data/testdata.csv", header=TRUE, stringsAsFactors = FALSE)
test_mirna_data %>% 
  mutate(mirna_class = factor("Mirtron")) %>%
  tibble::rownames_to_column() %>% 
  filter(!(as.numeric(rowname) %in% c(1,22,103,139,151,164,165,182,202))) %>% 
  mirna_features -> test_data

## PLOTS ####

#mirtronplots_mirna=mirnaplots(mirtron_mirna)
#canonicalplots_mirna=mirnaplots(canonical_mirna)
#testplots_mirna=mirnaplots(test_mirna)
mirna_means = mirna_data %>% 
  group_by(class) %>% 
  summarise_if(is.numeric, funs(mean))

plot_histogram = function (data, name, means){
  ggplot(data = data, aes_string(name, fill = "class", color = "class")) + geom_histogram(alpha = 0.3, bins = 30) +
  geom_vline(data = means, aes_string(xintercept=name, color="class"), linetype="dashed", size = 1.2) + 
  scale_fill_brewer(palette="Dark2") +
  scale_color_brewer(palette="Dark2") + 
  theme_minimal()+theme(legend.position="top")
}

plot_histogram = function (data, name, means){
  ggplot(data = data, aes_string(name, fill = "class", color = "class")) + geom_histogram(alpha = 0.3, bins = 30) +
    scale_fill_brewer(palette="Dark2") +
    scale_color_brewer(palette="Dark2") + 
    theme_minimal()+theme(legend.position="top")
}

#myplots <- lapply(mirna_means %>% select_if(is.numeric) %>% names, plot_histogram, data = mirna_data, means = mirna_means)
#multiplot(plotlist = myplots, cols = 5)
content = mirna_data %>% 
  select(mature5p_A, mature5p_C, mature5p_G, mature5p_U, class) %>% 
  rownames_to_column %>% gather("variable", "value", 2:5)

## STATISTICS ####
#source('mirnaplots.R')
P.values = data_frame()
sign=c('Not Significant','Significant')

tukey = function(x, classes) {
  sign=c('Not Significant','Significant')
  tukey_df=data_frame(input = x, classes = classes)
  means=tapply(tukey_df$input,tukey_df$classes, mean)
  ajuste <- lm( tukey_df$input ~ tukey_df$classes)
  tukey.alpha=0.01
  h=HSD.test(ajuste, 'tukey_df$classes',alpha=tukey.alpha)
  tukey = sign[length(levels(h$groups$groups))] #$M
  return(tukey)
}
classes = as.numeric(mirna_data$class == "Mirtron")
mirna_data %>% 
  summarise_if(is.numeric, funs("Wilcoxon - test" = wilcox.test(.[class == "Mirtron"], .[class == "Canonical"])$p.value,
                                "KS - test" = ks.test(.[class == "Mirtron"], .[class == "Canonical"])$p.value,
                                "Tukey - test" = tukey(., classes))) %>%
  gather %>%
  transmute(name = sub("\\_[^\\_]*$", "", key), test = sub(".*_", "", key), value = value) %>%
  spread(test, value) -> test_results

cat("\nStatistical tests\n")
print(test_results, n = Inf)

## Principal Component Analysis ####

ml_data = mirna_data %>% 
  dplyr::select(-c(contains("position"), contains("_U"), rowname, hairpin_name))
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

## CORRELATION ####
mir_cor = cor(ml_data %>% select(-class))
find_cor = findCorrelation(mir_cor, cutoff = 0.8, verbose = T, names = T)
ml_data = ml_data %>% select_if(!(names(.) %in% find_cor))
cor_df = as.data.frame(as.table(mir_cor)) %>%
  mutate(Freq = abs(Freq)) %>%
  arrange(desc(Freq)) %>%
  filter(Freq != 1) %>%
  filter(row_number() %%2 == 0) %>%
  rename(replace = c("Freq" = "Pearson"))

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
  mutate(Method = colnames(ml_data %>% select(-class))) %>%
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


g = ggplot(stepwise_results, aes(x = feature, F1)) + 
  geom_point() +
  scale_x_discrete(limits = stepwise_results$feature) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(g)

## VERIFICATION ####
no_best_features = which(stepwise_results$F1 == max(stepwise_results$F1))
final_data = ml_data %>% select(stepwise_results$feature[1:no_best_features], class)
x2=LogReg(final_data, folds)
pred=predict(x2$svm,ml_test_data) #predict on svm model
all = append(pred, ml_test_data$class)
conf_mat = table(as.character(pred), test_data$class)
accuracy = conf_mat[2] / sum(conf_mat)

## Print results ####
print(test_results, n = Inf)
print(x$results, n = Inf)
print(x2$results, n = Inf)
print(Singlef, n = Inf)
print(stepwise_results, n = Inf)
print(Boruta_results, n = Inf)
print(bind_cols(Boruta_results %>% filter(!grepl("shadow",Feature)), stepwise_results), n = Inf)
## Arrows ####
v1 = Boruta_results %>% filter(!grepl("shadow",Feature)) %>% select(Feature) %>% pull
v2 = stepwise_results %>% select(feature) %>% pull
o <- 0.07
DF <- data.frame(x = c(rep(1, length(v1)), rep(1.5, length(v2))),
                 x1 = c(rep(1 + o, length(v1)), rep(1.5 - o, length(v2))),
                 y = c(rev(seq_along(v1)), rev(seq_along(v2))),
                 g = c(v1, v2))
differences = diff(as.matrix(spread(DF,g,y)))[,-c(1:2)] %>% data.frame(change = .) %>% rownames_to_column() %>% rename(replace =  c("rowname" = "g"))
DF = join(DF, differences) %>% mutate(arr_color = ifelse(change > 0, "green", ifelse(change == 0, "yellow", "red")))
library(ggplot2)
library(grid)
ggplot(DF, aes(x=x, y=y, group=g, label = g)) +
  geom_path(aes(x=x1), arrow = arrow(length = unit(0.02,"npc")), 
            size=0.5, color=DF$arr_color) +
  geom_text(size=4) +
  theme_minimal() + theme(axis.title = element_blank(),
                          axis.text = element_blank(),
                          axis.ticks = element_blank(),
                          panel.grid = element_blank()) + xlim(c(0.9,1.6)) + 
  ggtitle("Boruta ranking vs SFS ranking") +
  theme(plot.title = element_text(hjust = 0.5))

