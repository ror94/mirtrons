mirnaplots = function(input){
attach(input)
plotlist=list()
hist(hairpin_length,breaks=40, col="red", xlim=c(50,100), main="Precusors length distribution",xlab="Length",ylab="Counts")
plotlist[[1]]=recordPlot()
hist(mature5p_length,breaks=10, col="red", xlim=c(16,27), main="5p mature length distribution",xlab="Length",ylab="Counts")
plotlist[[2]]=recordPlot()
hist(mature3p_length,breaks=10, col="red", xlim=c(16,27), main="3p mature length distribution",xlab="Length",ylab="Counts")
plotlist[[3]]=recordPlot()
slices1=c(mean(hairpin_A),mean(hairpin_C),mean(hairpin_G),mean(hairpin_U))
slices2=c(mean(mature5p_A),mean(mature5p_C),mean(mature5p_G),mean(mature5p_U))
slices3=c(mean(mature3p_A),mean(mature3p_C),mean(mature3p_G),mean(mature3p_U))
slices=data.frame(slices1,slices2,slices3)
pct <- round(slices)
lbls=c("A","C","G","U")

pie(slices$slices1, paste(paste(lbls,pct$slices1),"%",sep=""), main="Hairpin mean nucleotide composiotion")
plotlist[[4]]=recordPlot()
pie(slices$slices2, paste(paste(lbls,pct$slices2),"%",sep=""), main="5p mature mean nucleotide composiotion")
plotlist[[5]]=recordPlot()
pie(slices$slices3, paste(paste(lbls,pct$slices3),"%",sep=""), main="3p mature mean nucleotide composiotion")
plotlist[[6]]=recordPlot()
hist(Harpin_FE,breaks=20, col="red",xlim=c(-1,0),main="Free energy distribution",xlab="Free energy [kcal/mol]",ylab="Counts")
plotlist[[7]]=recordPlot()
hist(overhangnumber,breaks=20, col="red",xlim=c(-10,10),main="Overhang distribution",xlab="Overhang",ylab="Counts")
plotlist[[8]]=recordPlot()

return(plotlist)


}
