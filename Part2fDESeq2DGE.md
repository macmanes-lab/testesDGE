#Manuscript: Differential Gene Expression for TESTES


______

**Part 2) Differential Gene Expression Analysis**

#f. differential gene expression with tximport and DESeq2


**DEPENDENCIES:**

	source("https://bioconductor.org/biocLite.R")

	biocLite("tximport")

	library(tximport)

	install.packages("readr")

	library(readr)

	biocLite("DESeq2")

	library(DESeq2)


**import the gene ID matrix:**

	tx2gene <- read.delim("~/Desktop/NEWSALMONdge/NEWESTFinalMUS.csv", header=FALSE)
	head(tx2gene)


***SET DIRECTORY***

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


**Now read in files with tximport:**

	txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, read = read_tsv)
	head(txi.salmon$counts)


**Read in table for group names:**

	sampleTable <- read.csv("~/Desktop/mergedCounts/sampleTable.csv", header=TRUE, row.names=1)
	head(sampleTable)

*(Note: Dry = 2 ; Wet = 1)*

**Assign the condition column as a factor:**

	sampleTable$condition <- as.factor(sampleTable$condition) 
	

**Next, check if all true, and if all are true, then proceed:**

	rownames(sampleTable)==colnames(txi.salmon$counts)
	

**Prepare the entire matrix for analysis by condition:**

	dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~condition)


**Filter out the rows (transcript IDs) with less than 1 read:**

	dds <- dds[ rowSums(counts(dds)) > 1, ]

*(note that the count data is counts(dds))*


**Perform DEseq statistical analysis:**

	dds <- DESeq(dds)

**Results filtering based on normalized counts of each gene with a FDR value set at alpha= 0.05 for adjusted p value:** 

	res <- results(dds, alpha=0.05, pAdjustMethod = "BH")

	summary(res)


*RESULTS:*

out of 14216 with nonzero total read count

adjusted p-value < 0.05

LFC > 0 (up)     : 67, 0.47% 

LFC < 0 (down)   : 148, 1% 

outliers [1]     : 0, 0% 

low counts [2]   : 276, 1.9% 

(mean count < 1)

[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results


	total<-sum(res$padj < 0.05, na.rm=TRUE)

	total

[1] 215

*This means there are 215 significant results*


**Order the results from most to lease significant by adjusted p-value:**

	resOrdered <- res[order(res$padj),]


**All of the results were written to a file so that I could view them:**

	write.csv(resOrdered[1:total,],file="NEWresults.csv")


**Plot:**

	plotMA(resMLE, MLE=TRUE, main="unshrunken LFC", ylim=c(-4,4), colLine=NONE)
	abline(h = c(-2, 2), col = "blue")
	abline(h = c(-4, 4), col = "blue")
