####PANTHER analysis for evaluating median normalized cpm values between WET and DRY treatments for each gene (without significance cutoffs). This was executed on the edgeR DGE dataset in R:

	group.dif<-cpm(b)
	colnames(group.dif)<-group
	group.dif.1<-group.dif[,which(colnames(group.dif)=="1")]
	group.dif.2<-group.dif[,which(colnames(group.dif)=="2")]
	median.group.1<-apply(group.dif.1, 1, median)
	median.group.2<-apply(group.dif.2, 1, median)

	plot(median.group.1,median.group.2,ylim=c(0,5000),xlim=c(0,5000),ylab="Dry 	cpm",xlab="Wet cpm",main="Median")
	abline(0,1,col="red")
	
	correlation<-cor(medan.group.1,median.group.2)
	correlation
	correlation^2
	correlation.test<-cor.test(median.group.1,median.group.2)
	head(correlation.test) 

**The following code separates the relatively high DRY expressed genes (high.csv) from the relatively high WET expressed genes (low.csv) into two Ensembl ID gene lists.  After these two files were exported, I only retained the column in each file with the Ensembl IDs, and these two files were renamed highENSM.csv and lowENSM.csv, respectively, for upload into PANTHER.  These two files are available on DropBox.**

	high<-names(which((median.group.2/median.group.1)>1))
	write.csv(high,file="high.csv")
	low<-names(which((median.group.2/median.group.1)<1))
	write.csv(low,file="low.csv")

