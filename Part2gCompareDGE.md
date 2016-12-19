#Manuscript: Differential Gene Expression for TESTES


______

**Part 2) Differential Gene Expression Analysis**
	
	
#g. Comparison of Log2FC values for DGE in edgeR with DGE in DESeq2

***The prerequisite for this analysis is that the edgeR and DESeq2 analyses have already been completed in the current R session***
	

**Extract log2 FC values for each ensembleID in the et matrix from EdgeR:**

	output <- topTags(et, n=14218)

	write.table(output, file="diffexp_detags_edgeR.csv", sep = "," , row.names = TRUE)


*(Note: within et (my table), n = the total number of ensembl gene IDs)*


**Extract log2 FC values for each ensembleID in the res matrix from DESeq2:**

	write.table(res, file="diffexp_DESEQ2.csv", sep = "," , row.names = TRUE)

*(Note: my data table is res)*



**NEXT Extract data frame from the edgeR result object (from the list which is "output" select the .Data item in the list and take the first element of .Data, which is the results table, and convert this to a dataframe, named edgeRtable:**

	edgeRtable<-data.frame(output@.Data[1])


**Add rownames (Gene ID) as a column in each table:**

	edgeRtable$ID<-rownames(edgeRtable)

	res$ID<-rownames(res)


**Keep only the columns wanted and convert Fdeseq to dataframe:**

	Fdeseq<-as.data.frame(res[,c("log2FoldChange","ID")])

	FedgeR<-edgeRtable[,c("logFC","ID")]


**Merge final edgeR table with final DESeq2 table:**

	compare<-merge(Fdeseq,FedgeR,by="ID")


**Change the column names so that we label each lFC by the analysis:**

	names(compare)<-c("ID","D_l2FC","E_l2FC")


**linear regression of results (so that we can plot with a regression line):**

	regression<-lm(compare$E_l2FC~compare$D_l2FC)

	summary(regression)


Call:
	lm(formula = compare$E_l2FC ~ compare$D_l2FC)

Residuals:
     Min: -2.38746 
     1Q: -0.09001   
     Median: -0.00879       
     3Q: 0.09039      
     Max: 2.32470 
 

Coefficients:
               ***Estimate; Std. Error; t value; Pr(>|t|)***    
(Intercept)    0.305650   0.001669   183.1   <2e-16 ***

compare$D_l2FC 1.380585   0.008319   166.0   <2e-16 ***
	
	
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.1974 on 14214 degrees of freedom

Multiple R-squared:  0.6596,	Adjusted R-squared:  0.6596 

F-statistic: 2.754e+04 on 1 and 14214 DF,  p-value: < 2.2e-16



**Plot with regression line and axis labels, and regression line optional:**

	plot(compare$D_l2FC,compare$E_l2FC,ylab="EdgeR_lfc",xlab="DESeq2_lfc")
	abline(regression)


**Correlation of results:**

	correlation<-cor(compare$D_l2FC,compare$E_l2FC)

[1] 0.8121646
	
	correlation^2

[1] 0.6596113

	correlation.test<-cor.test(compare$D_l2FC,compare$E_l2FC)

Pearson's product-moment correlation

data:  compare$D_l2FC and compare$E_l2FC

t = 165.96, df = 14214, p-value < 2.2e-16

alternative hypothesis: true correlation is not equal to 0

95 percent confidence interval:

0.8064933 0.8176864

sample estimates:

cor 

0.8121646 
