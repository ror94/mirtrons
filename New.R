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
#install_github("vqv/ggbiplot")
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
library(ggrepel)
library(scales)
library(xtable)
#library(latticeExtra)
source("utils.R")
#source("mirna_features.R")
#source("overhangcount.R")
#source("loopscount.R")
#source("LogReg.R")
#source("BinEval.R")
#source("multiplot.R")


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

#rep(1.1*ymax, 2)

mirna_means = mirna_data %>% 
  group_by(class) %>% 
  summarise_if(is.numeric, funs(mean))

pos = cbind(mirna_means[,1],tbl_df(apply(mirna_means[,-1], 2, function(x) c(ifelse(x[1] > x[2], x[1] + 0.2*mean(x), x[1] - 0.2*mean(x)), 
                                       ifelse(x[2] > x[1], x[2] + 0.2*mean(x), x[2] - 0.2*mean(x))))))
pos = cbind(mirna_means[,1],tbl_df(apply(mirna_means[,-1], 2, function(x) c(ifelse(x[1] > x[2], 0, 1), 
                                                                            ifelse(x[2] > x[1], 0, 1)))))

plot_histogram = function (data, name, means, pos){
  g = ggplot(data = data, aes_string(name, fill = "class", color = "class")) + geom_histogram(alpha = 0.3, bins = 30) +
  geom_vline(data = means, aes_string(xintercept=name, color="class"), linetype="dashed", size = 1.2)
  ymax = max((ggplot_build(g))$data[[1]]$y)
  g + geom_label(data = means, aes_string(x = means %>% 
                                       select_(name) %>% 
                                       pull,
                                       hjust = pos %>% 
                                         select_(name) %>% 
                                         pull,
                                       size = 10
                                       #  = 0
                                       #color = class
                                     ), y = 1.05*ymax, vjust = 1, color = 'white',
                      size = 5, label = means %>% select_(name) %>% pull %>% format(digits = 3)) +
    scale_fill_manual(values = c( "#F8766D", "#619CFF")) +
    scale_color_manual(values = c( "#F8766D", "#619CFF")) +
    #scale_fill_brewer(palette="Dark2") +
    #scale_color_brewer(palette="Dark2") + 
  theme_minimal()+theme(legend.position="top")
}
plot_histogram(data = mirna_data, means = mirna_means, name = "mature5p_A", pos)

myplots <- lapply(mirna_means %>% 
                    select_if(is.numeric) %>% 
                    select(4,5,6,7,2,8,9,10,11,3,13,14,15,16,12,18,19,20,21,1,17,22,23,24,25) %>% 
                    names, tryCatch({plot_histogram}, error = function(e){print(name)}) , data = mirna_data, means = mirna_means, pos = pos)
multiplot(plotlist = myplots[1:25], cols = 5)
#content = mirna_data %>% 
#  select(mature5p_A, mature5p_C, mature5p_G, mature5p_U, class) %>% 
#  rownames_to_column %>% gather("variable", "value", 2:5)

#library(gridExtra)
#do.call("grid.arrange", c(myplots, ncol=5))

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
  summarise_if(is.numeric, funs("Wilcoxon test" = wilcox.test(.[class == "Mirtron"], .[class == "Canonical"])$p.value,
                                "KS - test" = ks.test(.[class == "Mirtron"], .[class == "Canonical"])$p.value,
                                "Tukey - test" = tukey(., classes))) %>%
  gather %>%
  transmute(name = sub("\\_[^\\_]*$", "", key), test = sub(".*_", "", key), value = value) %>%
  spread(test, value) %>% 
  inner_join(data.frame("Mirtron median" = apply(mirna_data 
                                               %>% filter(class == "Mirtron") 
                                               %>% select_if(is.numeric) , 2, FUN = median)) 
             %>% rownames_to_column(var = "name")) %>% 
  inner_join(data.frame("Canonical median" = apply(mirna_data 
                                                 %>% filter(class == "Canonical") %>% 
                                                   select_if(is.numeric) , 2, FUN = median)) %>% 
               rownames_to_column(var = "name")) %>%
  select(-2, -3) -> test_results
test_results$Mirtron.median = format(test_results$Mirtron.median, digits = 2 )
test_results$Canonical.median = format(test_results$Canonical.median, digits = 2 )
test_results$`Wilcoxon test`= format(as.numeric(test_results$`Wilcoxon test`), digits = 3 )
cat("\nStatistical tests\n")
print(test_results, n = Inf)


## Principal Component Analysis ####

ml_data = mirna_data %>% 
  dplyr::select(-c(contains("position"), contains("_U"), rowname, hairpin_name))
pca_data = dplyr::select(ml_data, -class)
pca=prcomp(pca_data, retx=TRUE, center=TRUE, scale=TRUE)
labels=factor(mirna_data$class)
grbiplot(pca, class = labels, scale_color_manual(values = c("#F8766D", "#619CFF")))


ml_test_data = test_data %>%
  dplyr::select(-c(contains("position"), contains("_U"), rowname, hairpin_name, interarm3p, interarm5p))
pca_test_data = dplyr::select(ml_test_data, -class)
pca_test_data = bind_rows(pca_data, pca_test_data)
pca2 = pca
pca2$x = scale(pca_test_data, pca$center, pca$scale) %*% pca$rotation
labels = factor(c(as.vector(mirna_data$class), rep('test',nrow(pca_test_data) - nrow(pca_data))))
grbiplot(pca2, class = labels, scale_color_manual(values = c("#F8766D", "#619CFF", "#00BA38")))


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
#mir_cor = cor(ml_data %>% select(-class))
#find_cor = findCorrelation(mir_cor, cutoff = 0.8, verbose = T, names = T)
#ml_data = ml_data %>% select_if(!(names(.) %in% find_cor))
#cor_df = as.data.frame(as.table(mir_cor)) %>%
#  mutate(Freq = abs(Freq)) %>%
#  arrange(desc(Freq)) %>%
#  filter(Freq != 1) %>%
#  filter(row_number() %%2 == 0) %>%
#  rename(replace = c("Freq" = "Pearson"))

## CLASSIFICATION ####
itnumber = 5
set.seed(27)
folds = generateCVRuns(ml_data$class, ntimes = 1, nfold = itnumber, stratified = TRUE)
cat("\n1x5-fold Cross-validation classification\n")
x=LogReg(ml_data, folds)
print(x$results)
pred=predict(x$svm, ml_test_data) #predict on svm model
ref = factor(c(rep("Mirtron",length(pred)), "Canonical"))[1:length(pred)]
cm = (confusionMatrix(pred,ref))$table
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
par(pty = "m")
b = plot(Bor, xaxt = "n", xlab = "", yaxt = "n")
axis(1, labels = FALSE, at = c(1:length(colnames(Bor$ImpHistory))))
# Plot x labs at default x position
text(x =  seq_along(colnames(Bor$ImpHistory)), y = par("usr")[3] - 2, srt = 45, adj = 1, cex = 0.7,
     labels = names(lz), xpd = TRUE)
axis(2, cex.axis=0.7, las = 1)


Boruta_results = data_frame(Feature = colnames(Bor$ImpHistory), Means = colMeans(Bor$ImpHistory)) %>%
  arrange(desc(Means))

## STEPWISE SVM ####
stepwise_results = stepwise(ml_data,folds)


g = ggplot(stepwise_results, aes(x = feature, F1)) + 
  geom_point() +
  scale_x_discrete(limits = stepwise_results$feature) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(aspect.ratio = 0.4)
print(g)

## VERIFICATION ####
no_best_features = which(stepwise_results$F1 == max(stepwise_results$F1))
final_data = ml_data %>% select(stepwise_results$feature[1:no_best_features], class)
x2=LogReg(final_data, folds)
pred2=predict(x2$svm, ml_test_data) #predict on svm model
cm2 = (confusionMatrix(pred2,ref))$table
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
print(cm)
print(cm2)

xtable(test_results ,digits=-5)
xtable(x$results ,digits = 3)
xtable(Singlef, digits = 3)
xtable(stepwise_results %>% select(-rowname), digits = 3)
xtable(Boruta_results, digits = 3)
xtable(x2$results ,digits= 3)
xtable(cm)
xtable(cm2)

## Arrows ####
v1 = Boruta_results %>% filter(!grepl("shadow",Feature)) %>% select(Feature) %>% pull
v2 = stepwise_results %>% select(feature) %>% pull
o <- 0.07
DF <- data.frame(x = c(rep(1, length(v1)), rep(1.5, length(v2))),
                 x1 = c(rep(1 + o, length(v1)), rep(1.5 - o, length(v2))),
                 y = c(rev(seq_along(v1)), rev(seq_along(v2))),
                 g = c(v1, v2))
differences = diff(as.matrix(spread(DF,g,y)))[,-c(1:2)] %>% data.frame(change = .) %>% rownames_to_column() %>% rename(g = rowname)
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

