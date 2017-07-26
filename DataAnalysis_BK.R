#The function performs robust Z-score normalization
#It subtracts the median and devides by the IQR
#Guide to Inteligent Data Analysis, p. 130
normalizeZRobust<-function(x,na.rm=T)
{
xNorm<-((x-median(x,na.rm=na.rm))/IQR(x,na.rm=na.rm))
}

#Function for robust univariate outlier detection 
#(Hampel, 1974) (Liu,et al. 2004) (Rousseeu, 1993)
outlierDet_Hampel<-function(x,c=1.4826, beta=3)
{
#c - consistency parameter - for normal distribution assumption 1.4826
#beta - number of standard deviations that define the outlier region
	med=median(x)
	MAD=c*median(abs(med-x))
	outlierRegion<-c(med-beta*MAD,med+beta*MAD)
	outliersID=which(x<outlierRegion[1] | x>outlierRegion[2])
	xClean<-x[-outliersID]
	list(xClean=xClean,outliersID=outliersID,outRange=outlierRegion)
}


#The function removes all rows with at least one NA
removeNA_Inf<-function(df)
{
	#PREPARE NA REMOVAL
	rm.na_rows=c()
	na.tab<-is.na(df)
	for (r in 1:dim(na.tab)[1]){
		if (!all(!na.tab[r,])){
			rm.na_rows=c(rm.na_rows,r)
		}
	}
	
	#PREPARE INF REMOVAL
	# rm.inf_rows=c()
	# inf.tab<-is.infinite(as.matrix(df))
	# for (r in 1:dim(inf.tab)[1]){
	#   if (!all(!inf.tab[r,])){
	#     rm.inf_rows=c(rm.inf_rows,r)
	#   }
	# }
	rm.inf_rows<-unlist(lapply(df,function(x){which(is.infinite(x))}))
	rm.rows=unique(c(rm.inf_rows,rm.na_rows))
  if (length(rm.rows))
  {
    return (df<-list("data"=df[-rm.rows,],'rmInd'=rm.rows, "naInd"=rm.na_rows,"infInd"=rm.inf_rows))
  } else
    df<-list("data"=df,"rmInd"=NULL, 'naInd'=NULL, 'infInd'=NULL)
}

#The function checks the correlations between features in a data table
#It lists all that are greater than a given input threshold, it skips cor=1
BK_checkCorr<-function(data, threshold=0.8)
{
  #CHECK DATATABLE FOR HIGH CORRELATIONS
  corMat<-cor(data)
  corInd=which(corMat>0.8 & corMat !=1,arr.ind = T)
  highCorr=data.frame()
  if (dim(corInd)[1]>=1)
  {
    for(i in 1:dim(corInd)[1]) 
    {
      print(c(names(data)[corInd[i,1]],names(data)[corInd[i,2]],
            corMat[corInd[i,1],corInd[i,2]],
            corInd[i,1],corInd[i,2]))
    }
  }
  corInd
}




#The function constructs a data frame based on given formula and input data
#It moves the response variable to the first column
BK.expand.formula<-function(formula,data)
{
  f<-as.formula(formula)
  data <- model.frame(data)
  y<-data.frame(data[,all.vars(f)[1]])
  names(y)<-all.vars(f)[1]
  row.names(y)<-row.names(data)
  x<-data.frame(model.matrix(f,data))
  res=data.frame(y,x[,-1])
}

#The function generates histgrams for all atrributes in dataframe
#Both normal values and their natural logarythims are ploted
BK_plotHists<-function(df,width=1240,height=720)
{
	namesList<-names(df)
	tmpDir<-paste("tmp_",round(abs(100*rnorm(1))),sep="")
	dir.create(tmpDir)
	
	for (i in seq(df))
	{			
		filePath<-paste(tmpDir,"/hist_",namesList[i],".png",sep="")	
		png(file=filePath,width=1240,height=720)
		par(mfcol=c(1,2))
		hist(df[[i]],xlab=namesList[i], main=namesList[i])
		hist(log(df[[i]]),xlab=namesList[i],main=paste("ln(",namesList[i],")",sep=""))
		dev.off()
	}
	
}

#The function adds a natural logarithm transformation of parametr X to source data
addLog<-function(df,attrNames)
{
	dfNames<-names(df)
	for (i in attrNames){
		if(!is.na(pmatch(i,dfNames)))	
		{
			logTrans<-data.frame(log(df[[i]]))
			names(logTrans)<-paste("log",i,sep="")
			df<-data.frame(df,logTrans)
		}
	}
df

}



BK_PCA_bplot_sdev<-function(input,center=TRUE,scale.=FALSE)
{
#The function performs PCA analysis and visualises the results using 
#biplot and a barplot with component standard deviations

	res<-prcomp(input,center=center)
	dev.new()
	par(mfcol=c(1,2))
	biplot(res,main="Component 1 and 2, original atributes directions")
	plot(res, main="Component standard deviations")
	res
}
BK_clv<-function(data,clust,dist='euclidean',ValDistIntra='centroid',ValDistInter='centroid')
{
#The function performs validation of clusters formed in the cluster - cluster asignment matrix. 
#Each column in cluster should be a vector of cluster labels. Varaiable "data" is the dataframBKe with object attributes.
#The function returns Davies-Bouldain and Dunn indices along with cls.scatt.data containing various intra
#and inter-cluster distances
	require('clv')

	f1<-function(x,y,z)
	{
		cls.scatt.data(clust=x, data=y, dist=z)
	}

	clDistData<-lapply(data.frame(clust),f1,y=data,z=dist)
	
	Davies.Bouldin<-sapply(clDistData,clv.Davies.Bouldin,ValDistIntra,ValDistInter)
	Dunn<-sapply(clDistData,clv.Dunn,ValDistIntra,ValDistInter)

	list("clDistData"=clDistData,"Davies.Bouldin"=Davies.Bouldin,"Dunn"=Dunn)
}

byGroup_concat<-function(Group_h2h)
{
	listed_p.values<-matrix(0,ncol=3)
	atr_names=names(Group_h2h)
	for (atr in 1:length(Group_h2h))
	{
		mat=as.matrix(Group_h2h[[atr]])
		row.names(mat)<-rep.int(atr_names[atr],dim(mat)[1])
		listed_p.values<-rbind(listed_p.values,mat)
				
	}
	storage.mode(listed_p.values)<-"double"
	listed_p.values=listed_p.values[-1,]
	listed_p.values=listed_p.values[order(listed_p.values[,"p.value"]),]
}

BK_head2head<-function(data,clust,p.adjust="bonferroni",alpha=0.01)
{
#The function first performs a Kruskal-Walis rank group test. If there are significant differences,
#it performs performs pairwise Mann-Whitney tests.

KruskalW_res<-lapply(data,kruskal.test,clust)#Kruskal-Walis test
KruskalW_p.values<-data.frame(lapply(KruskalW_res,f<-function(x){x$p.value}))

MannW_res<-lapply(data,pairwise.wilcox.test,clust,p.adjust=p.adjust)#Mann-Whitney or Wilcox pairwise rank sum test
MW_rm<-KruskalW_p.values>alpha


MannW_res[unlist(MW_rm)]=NULL
MannW_mat=lapply(MannW_res,f<-function(x){x$p.value})

MannW_list<-lapply(MannW_mat,head2head_extract_p.values)

byGroup_h2h<-byGroup_extract_p.values(MannW_list,unique(clust))

byGroup_listed<-lapply(byGroup_h2h,byGroup_concat)

list(KruskalW_p.values=KruskalW_p.values,MannW_p.values=MannW_list, byGroup_h2h=byGroup_h2h,byGroup_listed=byGroup_listed)
}

head2head_extract_p.values<-function(h2hMat)
{
 # print(h2hMat)
	p.value<-data.frame(h2hMat)#[,1]
	p.value<-p.value[!is.na(p.value)]
	id_names<-t(combn(unique(c(colnames(h2hMat),row.names(h2hMat))),2))
  #print(id_names)
 # print(p.value)
	res<-data.frame(p.value,id_names)
	
}

byGroup_extract_p.values<-function(p.values_list,groups)
{
	res=list()
	for (gr in groups)
	{
		res[paste(c("gr",gr),collapse="")]<-list(lapply(p.values_list,extractGroup,gr))
	}
	res
}

extractGroup<-function(p.values,gr)
{
	gr_p.values<-p.values[(p.values$X1==gr | p.values$X2==gr),]
	gr_p.values<-gr_p.values[order(gr_p.values$p.value),]
}
BK_binGenData<-function(genData)
{
	genData.bin<-data.frame(lapply(genData,f<-function(x){
		m=matrix(0,nrow=length(x),ncol=length(levels(x)))
		for (i in 1:length(x)) 
			{
				m[i,x[i]]=1
		}
		m
	}))
	row.names(genData.bin)<-row.names(genData)
	genData.bin
}

BK_lm<-function(data,formula=NULL,K=1,crossVal.ID=NULL,method="lm",init=NULL)
{
	require(MASS)
	require(sfsmisc)
	
  
	#checking formula
	if (is.null(formula))
	{
		lm.formula<-formula(data)
		print(lm.formula)
	}
  #creating the response vector
  response.name<-as.character(formula)[2]
	response=data[,response.name]
  response.no=which(names(data)==response.name)
	#preparing K-fold crossValidation
	if (is.null(crossVal.ID))
	{
		crossVal.ID<-BK_sample_crossVal(dim(data)[1],K)
	}
	else K=length(crossVal.ID)
	
	#buildnig models
	model.list=list()
	if (method=="rlm" & is.null(init)) 
	{
		print("Initializing coefficients with OLS")
		init<-lm(formula,data=data[unlist(crossVal.ID[-1]),])
		init<-coefficients(init)
	}
	
	
	for (i in 1:K )
	{
		cv.subset<-unlist(crossVal.ID[-i])
		if (method=="lm") model.list[[i]]<-lm(formula,data=data[cv.subset,])
		else if(method=="rlm") model.list[[i]]<-rlm(formula,data=data[cv.subset,],init=init,maxit=30)
			
	}
	names(model.list)<-paste("m",1:K,sep="")

	#Model evaluation

	if (method!="lm")
	{
		print("Robust Linear Model")
		#p.values for robust lm
		p.list<-lapply(model.list,BK_single_robust_ftest)
		anova.table<-NULL
		a.sig.pass<-NULL
		
	} else
	{
		print("Ordinary Linear Model")
		#p.values for standard lm
		p.list<-lapply(model.list,f<-function(x){summary(x)$coefficients[,4]}) #pobranie p.valu dla t-testu
		anova.table<-lapply(model.list,f<-function(x){anova(x)[-dim(data)[2],5]}) # pobranie p.value dla F-testu
		anova.table<-matrix(unlist(anova.table),ncol=K)
		row.names(anova.table)<-colnames(data[,-1])
		#print(dim(anova.table))
		#print(dim(data))
	}
  
  RSS=matrix(0,nrow=K)
  MSE=matrix(0,nrow=K)
  R2=matrix(0,nrow=K)
  for (i in 1:K)
  {
    cv.test.ID<-unlist(crossVal.ID[i])
	  k.pred <- predict(model.list[[i]],newdata= data[cv.test.ID,])
	  RSS[i] <- sum((k.pred-response[cv.test.ID])^2)
	  MSE[i] <- RSS[i]/length(response[cv.test.ID])
	
	  TSS <- sum((response[cv.test.ID] - mean(response[cv.test.ID]))^2)
	  R2[i] <- 1 - RSS[i]/TSS
  }
  cv.MSE<-sum(MSE)/K # calculating crossvalidation test MSE
  cv.R2<-sum(R2)/K   # calculating crossvalidation test R2
  
  cv.test<-list(cv.MSE=cv.MSE,cv.R2=cv.R2,MSE=MSE,R2=R2)
  
  
  
	p.table<-matrix(unlist(p.list),ncol=K)
	print(dim(p.table))
	print(dim(data[,-1]))
	row.names(p.table)<-c("(Intercept)", colnames(data[,-1]))
	
	
	#Calculating significnce test pass
	sig.pass=data.frame(1:dim(p.table)[1])
	if(!is.null(anova.table))a.sig.pass=data.frame(1:dim(anova.table)[1])
	alpha=c(0.05, 0.01, 0.001, 0.0001)
	for (a in alpha)
	{
		sig.pass=data.frame(sig.pass,rowSums(p.table<a))
		if (!is.null(anova.table))a.sig.pass=data.frame(a.sig.pass,rowSums(anova.table<a))
	}
	sig.pass<-sig.pass[-1]
	row.names(sig.pass)<-c("(Intercept)", colnames(data[,-1]))
	colnames(sig.pass)<-paste("<", alpha, sep="")
	if(!is.null(anova.table))
	{
		a.sig.pass<-a.sig.pass[-1]
		row.names(a.sig.pass)<-c(colnames(data[,-1]))
		colnames(a.sig.pass)<-paste("<", alpha, sep="")
	}
	
	
	#Extracting coefficients
	coef.table<-lapply(model.list,f<-function(x){summary(x)$coefficients[,1]})
	coef.table<-matrix(unlist(coef.table),ncol=K)
	row.names(coef.table)<-c("(Intercept)", colnames(data[,-1]))
  
  #Calculating mean coeff
  coef.mean<-data.frame(apply(coef.table,1,mean))
  coef.sd<-data.frame(apply(coef.table,1,sd))
  coef.ttest<-apply(coef.table,1,t.test)
	coef.ttest<-t(data.frame(lapply(coef.ttest,function(x){return(c(x$p.value,x$conf.int[1:2]))})))
  coef.mean<-data.frame(coef.mean,coef.sd,coef.ttest)
  names(coef.mean)<-c("Mean.coef","SD","p.value","CI_low","CI_high")
  format(coef.mean,digits=4)
  

	
	list(model.list=model.list,p.values=p.table, sig.pass=sig.pass, anova.pv=anova.table, anova.sig.pass=a.sig.pass, coef.mean=coef.mean, coef=coef.table, cv.test=cv.test, crossVal.ID=crossVal.ID)
}
BK_sample_crossVal<-function(N,K)
{
  #The function returns a list of index sets for a crossvalidation study. 
  #It takes as input parameters: N - number of instances 
  #and K - number of crossvalidation sets
  #EXAMPLE
  #crossVal.ID<-BK_sample_crossVal(dim(data)[1],K)
	crossVal.ID=list()
	n=round(N/K)
	shuffled.ID<-sample(1:N,N)
	if (K>1)
	{
		for (k in 1:(K-1))
		{
			crossVal.ID[[k]]=shuffled.ID[((k-1)*n+1):(k*n)]
		}
	}
	crossVal.ID[[K]]<-shuffled.ID[((K-1)*n+1):N]
	names(crossVal.ID)<-paste("cv.",c(1:K),sep="")
	crossVal.ID
}
BK_single_robust_ftest<-function(lm.model)
{
#Calculates the F-test p.values using f.robftest function from sfsmisc package
	robustF=list()
	for (i in names(coefficients(lm.model)))
	{
		robustF[i]=f.robftest(lm.model,var=i)$p.value
	}
	robustF
}

BK_cv_step<-function(K=10,formula,data,scope,direction="backward")
{
  N=dim(data)[1]
  p=dim(data)[2]-1
  crossVal.ID<-BK_sample_crossVal(N=N,K=K)
  model.list=list()
  for (i in 1:K )
  {
    cv.subset<-unlist(crossVal.ID[-i])
    print(length(cv.subset))
    cv.data=data[cv.subset,]
    model.list[[i]]<-step(lm(formula=formula,data=cv.data),scope=scope,direction=direction)
    #model.list[[i]]<-step(lm(formula=formula,data=data[cv.subset,]),scope=scope,direction=direction)
  }
  
  cv_eval<-BK_eval_cv_step(model.list=model.list,data=data,crossVal.ID=crossVal.ID,K=K)
  
}
BK_eval_cv_step<-function(model.list,data,crossVal.ID,K=K)
{
  #p.values for standard lm
  p=dim(data)[2] #number of predictors
  
  p.model<-lapply(model.list,f<-function(x){summary(x)$coefficients[,4]}) #extracting p.value of the t-test
  p.full=matrix(1,ncol=K,nrow=p)
  rownames(p.full)=c("(Intercept)",colnames(data[,-1]))
  colnames(p.full)=paste("m",1:K,sep="")
  print(dim(p.full))
  for (i in 1:K)
  {
    for (j in 1:length(p.model[[i]]))
   {
      name.model<-names(p.model[[i]])[j]
      name.full<-gsub(':','.',name.model)
      print(name.model)
      print(name.full)
      #print(p.table[name,i])
      name.full=tryCatch({p.full[name.full,i]
                      name.full},
                              error=function (cond) {
                              name.full=paste(strsplit(name.full,'[.]')[[1]][c(2,1)],collapse='.')
                              #p.table[name_p.table,i]=p.list[[i]][name]
                              return(name.full)})
      print(name.full)
      print(name.model)
      p.full[name.full,i]=p.model[[i]][name.model]
      
    }
  }
  print(p.full)
   #Calculating significnce test pass
  sig.pass=data.frame(1:dim(p.full)[1])
  alpha=c(1,0.05, 0.01, 0.001, 0.0001)
  for (a in alpha)
  {
    sig.pass=data.frame(sig.pass,rowSums(p.full<a))
  }
  sig.pass<-sig.pass[-1]
  row.names(sig.pass)<-c("(Intercept)", colnames(data[,-1]))
  colnames(sig.pass)<-paste("<", alpha, sep="")

  list(model.list=model.list,p.values=p.full, sig.pass=sig.pass, crossVal.ID=crossVal.ID)
# 
#   anova.list<-lapply(model.list,f<-function(x){anova(x)[-dim(data)[2],5]}) # extracting p.value of the F-test
#   
#   
#   anova.table<-matrix(unlist(anova.table),ncol=K)
#   row.names(anova.table)<-colnames(data[,-1])
#   #print(dim(anova.table))
#   #print(dim(data))
# 
#   p.table<-matrix(unlist(p.list),ncol=K)
#   print(dim(p.table))
#   print(dim(data[,-1]))
#   row.names(p.table)<-c("(Intercept)", colnames(data[,-1]))
# 
# 
#   #Calculating significnce test pass
#   sig.pass=data.frame(1:dim(p.table)[1])
#   if(!is.null(anova.table))a.sig.pass=data.frame(1:dim(anova.table)[1])
#   alpha=c(0.05, 0.01, 0.001, 0.0001)
#   for (a in alpha)
#   {
#     sig.pass=data.frame(sig.pass,rowSums(p.table<a))
#     if (!is.null(anova.table))a.sig.pass=data.frame(a.sig.pass,rowSums(anova.table<a))
#   }
#   sig.pass<-sig.pass[-1]
#   row.names(sig.pass)<-c("(Intercept)", colnames(data[,-1]))
#   colnames(sig.pass)<-paste("<", alpha, sep="")
#   if(!is.null(anova.table))
#   {
#     a.sig.pass<-a.sig.pass[-1]
#     row.names(a.sig.pass)<-c(colnames(data[,-1]))
#     colnames(a.sig.pass)<-paste("<", alpha, sep="")
#   }
# 
# 
#   #Extracting coefficients
#   coef.table<-lapply(model.list,f<-function(x){summary(x)$coefficients[,1]})
#   coef.table<-matrix(unlist(coef.table),ncol=K)
#   row.names(coef.table)<-c("(Intercept)", colnames(data[,-1]))
# 
#   #Calculating mean coeff
#   coef.mean<-data.frame(apply(coef.table,1,mean))
#   coef.sd<-data.frame(apply(coef.table,1,sd))
#   coef.ttest<-apply(coef.table,1,t.test)
#   coef.ttest<-t(data.frame(lapply(coef.ttest,function(x){return(c(x$p.value,x$conf.int[1:2]))})))
#   coef.mean<-data.frame(coef.mean,coef.sd,coef.ttest)
#   names(coef.mean)<-c("Mean.coef","SD","p.value","CI_low","CI_high")
#   format(coef.mean,digits=4)
#   
#   #Calculating test  MSE and test R2
#   
#   list(model.list=model.list,p.values,p.table, sig.pass=sig.pass, anova.pv=anova.table, anova.sig.pass=a.sig.pass, coef.mean=coef.mean, coef=coef.table, crossVal.ID=crossVal.ID)
}

BK_gg2_biplot <- function(PC, x="PC1", y="PC2",gr.colour=1,gr.shape=1,label=NULL,size=6,fontface = "bold",limits=NULL,attr_text_size=3) {
  
  #Addapted from: http://stackoverflow.com/questions/6578355/plotting-pca-biplot-with-ggplot2
  # PC being a prcomp object
  
  if(is.null(label)) {
    data <- data.frame(obsnames=row.names(PC$x), PC$x,gr.colour,gr.shape)
    plot <- ggplot(data, aes_string(x=x, y=y)) + geom_point(alpha=1, size=size, aes(shape=gr.shape,colour=gr.colour))
    plot <- plot + geom_hline(aes(yintercept=0), size=.2) + geom_vline(aes(xintercept=0), size=.2)
    }
  else{
    data <- data.frame(obsnames=row.names(PC$x), PC$x,gr.colour,gr.shape,label)
    plot <- ggplot(data, aes_string(x=x, y=y)) + geom_text(alpha=1, size=size, fontface=fontface, aes(colour=gr.colour,label=label))
    plot <- plot + geom_hline(aes(yintercept=0), size=.2) + geom_vline(aes(xintercept=0), size=.2)
  }
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = attr_text_size, vjust=1, color="red")
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")
  if (!is.null(limits)){
    plot<-plot+expand_limits(x=limits$x,y=limits$y)
  }

}

# Multiple plot function
#Copied from: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
BK_multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL,labels=NULL,fsize=40,fcol='violet') {
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
      if (!is.null(labels)){
        grid.text(labels[i],x=0.15,y=0.9,vp=viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col),gp=gpar(fontsize=fsize, col=fcol))
      }
    }
  }
}

#The function calculates basic measures for assesing regresion (residuals, rse, R2, adjusted R2)
#It takes a vector of targets, a data.frame of predictions for assessment and optionally p - 
#number of explonatory variables for calulating adjusted R2
BK.assessRegression<-function(targets,pred_df,p=NULL)
{
  #INIT adjustedR2 in case p - number of predictors is not supplied
  adjustedR2=NULL
  n=length(targets)
  #RESIDUALS
  e=targets-pred_df
  eSummary=summary(e)
  #RESIDUAL STANDARD ERROR
  rse=sqrt(colSums(e^2)/(length(targets)-2))
  
  RSS<-colSums((e)^2)
  
  TSS <- sum((targets - mean(targets))^2)
  R2 <- 1 - RSS/TSS
  if (!is.null(p)){
    adjustedR2<-R2-(1-R2)*(p-1)/(n-p)
  }
  #Correlation coeficient
  #r2=apply(pred_df,2,function(x,y,z){cor(x,y,method=z)},targets,cor_method)
  
  list(rse=rse,R2=R2,adjustedR2=adjustedR2,e=e,eSummary=eSummary)
}