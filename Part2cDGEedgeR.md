#Manuscript: Differential Gene Expression for TESTES


______

**Part 2) Differential Gene Expression Analysis**
	
#c. Differential Gene Expression (DGE) with tximport and edgeR:


**R DEPENDENCIES:**

	source("https://bioconductor.org/biocLite.R")

	biocLite("tximport")

	library(tximport)

	install.packages("readr")

	library(readr)

	biocLite("edgeR")

	library(edgeR)


**import the gene ID matrix:**

	tx2gene <- read.delim("~/Desktop/NEWSALMONdge/NEWESTFinalMUS.csv", header=FALSE)	
	head(tx2gene)


**SET DIRECTORY:**

	dir <- setwd("/Users/lauren/Desktop/NEWSALMONdge/test")
	getwd()
	file.path(dir, "salmon", "quant.sf")
	dir <- setwd("/Users/lauren/Desktop/NEWSALMONdge/test")
	getwd()
	file.path(dir, "salmon", "quant.sf")

*This is what should be returned by these commands:*
      "/Users/lauren/Desktop/NEWSALMONdge/test/salmon/quant.sf"  
      

**SET samples in downloaded order:**


	samp <- c("salmon_13T", "salmon_102T", "salmon_209T", "salmon_265T", "salmon_343T", 	"salmon_344T", "salmon_349T", "salmon_355T", "salmon_366T", "salmon_376T", 	"salmon_381T", "salmon_382T", "salmon_383T", "salmon_384T", "salmon_400T", 	"salmon_888T","salmon_999T", "salmon_1357T", "salmon_1358T", "salmon_1359T", 	"salmon_2322T", "salmon_3333T") 
	samp

	file.path(dir, "salmon", samp, "quant.sf")
	files <- file.path(dir, "salmon", samp, "quant.sf")
	files

	names(files) <- samp
	names(files)


**NOW read in files with tximport:**

	txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, read = read_tsv)
	head(txi.salmon$counts)


**Quantification with edgeR:**

	cts <- txi.salmon$counts
	normMat <- txi.salmon$length
	normMat <- normMat/exp(rowMeans(log(normMat)))


	o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
	y <- DGEList(cts)
	y$offset <- t(t(log(normMat)) + o)


**Assign treatment groups:**

	group <- c(2,1,2,2,2,1,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1)
	group


**Make DGEList:**

	y <- DGEList(counts=y, group=group)
	y


**Scaling plot:**

	plotMDS.DGEList(y , main = "MDS Plot for Count Data", labels = colnames(y))


**Estimate Common Dispersion:**

	cds <- estimateCommonDisp(y)
	names(cds)

	cds$common.dispersion

*0.5330103  (this is the value returned for common.dispersion)*


**Estimate Tagwise Dispersions:**

*n is the integer value should be the nearest integer to this eqn: 50/(#samples - #groups) = 50/(22-2) = 50/20 ~ 2*

	cds <- estimateTagwiseDisp(cds, prior.n=2)
	names(cds)
	summary(cds$tagwise.dispersion)

*Values returned:* 

Min: 0.1362 
	
1st Qu: 0.4216 
 
Median: 0.5196
 
Mean: 0.5302 

3rd Qu: 0.6228 
 
Max. 2.8920  
 

**I tested whether doing the common & tagwise dispersions a different way that does not incorporate calculating the tagwise dispersion effected the results, and it did not:**


	b <- estimateCommonDisp(y)
	b <- estimateTagwiseDisp(b)
	et <- exactTest(b)
	topTags(et)

	detags <- rownames(topTags(et, n=20))
	cpm(b)[detags,]
	summary(de <- decideTestsDGE(et, p=0.05, adjust="BH"))
	detags <- rownames(b)[as.logical(de)]



**Mean Variance Plot:**

	meanVarPlot <- plotMeanVar(cds , show.raw.vars=TRUE,
                           	show.tagwise.vars=TRUE,
                           	show.binned.common.disp.vars=FALSE,
                           	show.ave.raw.vars=FALSE,
                           	dispersion.method="qcml" , NBline=TRUE,
                           	nbins=100,
                           	pch=16,
                           	xlab="Mean Expresion (Log10 Scale)",
                           	ylab="Variance (Log10 Scale",
                           	main="Mean-Variance Plot")


**Exact Test simplified version:**

	et <- exactTest(cds, pair=c("1","2"))
	summary(de <- decideTestsDGE(et, p=0.05, adjust="BH"))
	 

*RESULTS:*

   [,1]
    
-1     7

0  14203

1      8


*this means that 7 are expressed more highly in wet (-1)*
 
*this means that 8 are more highly expressed in dry (1)*

*Thus, there are 15 total Statistically significant genes*



**Because there are 15 significant results, just ask for all of them:**

	topTags(et, n=15)


*When my fold change is positive it is more highly expressed in dry (dry=2, wet=1)*

*Results*

Comparison of groups:  2-1 

   
   logFC ;    logCPM  ;    PValue ;   FDR
                      
ENSMUSG00000079019.2  -4.353530  1.6501524 4.093742e-13 5.820482e-09
	
ENSMUSG00000001768.15  3.085565  1.0062911 2.054417e-11 1.460485e-07

ENSMUSG00000054200.6  -3.733526  0.6187631 3.830859e-10 1.815572e-06

ENSMUSG00000025479.9   2.970609  3.0005331 2.242635e-08 7.971446e-05

ENSMUSG00000026435.15 -2.447967  2.4467916 4.416182e-07 1.130067e-03

ENSMUSG00000025020.11 -2.231147  1.7702447 5.733518e-07 1.130067e-03

ENSMUSG00000020427.11  2.681363  3.8866441 6.311282e-07 1.130067e-03

ENSMUSG00000019997.11  2.314055  3.2347652 6.807368e-07 1.130067e-03

ENSMUSG00000031170.14 -2.420675  2.5776205 7.153326e-07 1.130067e-03

ENSMUSG00000040170.13  1.950605  0.7526918 1.212177e-06 1.723474e-03

ENSMUSG00000023915.4   1.534036  1.2898632 1.565007e-05 2.022843e-02

ENSMUSG00000052974.8   2.077249  0.6465702 1.908666e-05 2.261451e-02

ENSMUSG00000030830.18 -2.179880  1.6659841 3.081597e-05 3.370319e-02

ENSMUSG00000027901.12  2.491578 -0.6204254 4.704997e-05 4.778261e-02

ENSMUSG00000032554.15 -2.066314  3.2865685 5.111751e-05 4.845258e-02


*I looked up the ensembl IDs for these significant results*



**plot:**

	detags <- rownames(cds)[as.logical(de)]
	plotSmear(cds, de.tags=detags, ylim=c(-6,6), xlim=c(0,20), frame.plot="false")
	abline(h=2, col='blue')
	abline(h=-2, col='blue')

--

***Analysis by BoxPlot of CPM values for differences by TRT in nine genes of interest***



ENSMUSG00000079019.2  Insl3   Insulin-like 3
	
ENSMUSG00000001768.15   Rin2    Ras and Rab Interactor 2
	
ENSMUSG00000054200.6  Ffar4   Free Fatty Acid Receptor 4
	
ENSMUSG00000026435.15  Slc45a3 solute carrier family 45, member 3
	
ENSMUSG00000020427.11   Igfbp3  Insulin-like growth factor binding protein 3
	
ENSMUSG00000019997.11    Ctgf    Connective Tissue growth factor
	
ENSMUSG00000031170.14  Slc38a5 solute carrier family 38, member 5
	
ENSMUSG00000030830.18   Itgal   integrin alpha L
	
ENSMUSG00000032554.15   Trf     Transferrin


**I graphed all nine of these genes (their cpm values for WET vs. DRY):**

	tag.name<-c("Rin2", "Igfbp3", "Ctgf", "Insl3", "Ffar4", "Slc45a3", "Slc38a5", "Itgal", "Trf")
	tag.keep<-c("ENSMUSG00000001768.15", "ENSMUSG00000020427.11", "ENSMUSG00000019997.11",
			"ENSMUSG00000079019.2", "ENSMUSG00000054200.6", "ENSMUSG00000026435.15",
            "ENSMUSG00000031170.14","ENSMUSG00000030830.18","ENSMUSG00000032554.15")
	group.dif<-cpm(b)
	group.dif<-group.dif[tag.keep,]
	colnames(group.dif)<-group
	group.dif<-melt(group.dif)
	names(group.dif)<-c("ID","group","cpm")
	group.dif$group[group.dif$group==1]<-"WET"
	group.dif$group[group.dif$group==2]<-"DRY"
	par(mfrow=c(3,3))
	for(i in 1:length(tag.keep)){
  		boxplot(cpm~group,group.dif[group.dif$ID==tag.keep[i],],
          main=tag.name[i],ylab="count (cpm)",frame.plot=F,col=(c("brown","light blue")))
	}

