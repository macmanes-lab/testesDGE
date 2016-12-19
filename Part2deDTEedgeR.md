#Manuscript: Differential Gene Expression for TESTES


______

**Part 2) Differential Gene Expression Analysis**

#d. MERGE QUANT FILES
	
***The purpose of this merging step is so that the data is compiled properly for DTE***
	
*All of the quant files are in one folder, and I used the following MacManes script to merge the quant files into one file.*

The merge_quants.py script to execute this merging step was saved in the salmon folder, where all of the quant outputs are.


**The following is the merge_quants.py script:**

	#!/usr/bin/env python
	#USAGE: python merge_quants.py dry/ dry.counts


	import os
	import sys
	from os import path
	import pandas as pd



	listofsamples = os.listdir(sys.argv[1])
	quants = {}
	data = None
	for sample in listofsamples:
    	if os.path.isdir(sys.argv[1]+sample):
        	if os.path.isfile(sys.argv[1]+sample+"/quant.sf"):
            	quant_file = sys.argv[1]+sample+"/quant.sf"
            	data=pd.DataFrame.from_csv(quant_file,sep='\t')
            	numreads = data['NumReads']
            	quants[sample] = numreads

	counts = pd.DataFrame.from_dict(quants)
	counts.set_index(data.index,inplace=True)
	counts.to_csv(sys.argv[2]



**The command to run the script was:**

	python merge_quants.py /mnt/data3/lauren/TESTESdge/NEWmaps/salmon/ NEWmerged.counts

*This file was also downloaded and named: NEWmerged_counts.csv*


The data appear with the transcript ID as the left most column, and then each sample is an additional column, with the raw transcript counts within each sample for that transcript ID.  

Thus, each row contains the unique transcript ID with the counts of transcripts (raw counts, NOT TPMs) for each sample.  



#e. differential transcript expression with edgeR


**DEPENDENCIES:**

	source("https://bioconductor.org/biocLite.R")

	biocLite("tximport")

	library(tximport)

	install.packages("readr")

	library(readr)

	biocLite("edgeR")

	library(edgeR)


**My dataset is NEWmerged.counts (and it is csv form) and I import it below:**


	x <- read.csv("~/Desktop/NEWSALMONdge/NEWmerged.counts", header=TRUE, row.names=1)
	head(x)


**Because only 1 sample needs a cpm >1 to keep the transcript ID , this is logical to filter rows by:**

	z <- x
	keep <- rowSums(cpm(z)>1) >= 1
	y <- z[keep,]
	dim(y)
	y[y==0]<-1


**Assign groups:**

	group <- factor(c(1,1,1,1,2,2,2,2,1,2,1,1,1,1,2,1,2,2,2,1,2,2))

**DGE Analysis:**

	b <- DGEList(counts=y,group=group)
	b <- calcNormFactors(b)
	b <- estimateCommonDisp(b)
	b <- estimateTagwiseDisp(b)
	et <- exactTest(b)
	topTags(et)

	detags <- rownames(topTags(et, n=20))
	cpm(b)[detags,]
	summary(de <- decideTestsDGE(et, p=0.05, adjust="BH"))
	detags <- rownames(b)[as.logical(de)]


*These are my Results:*

   [,1]
    
-1    45

0  69021

1     21


*45 higher in WET (-1)*
*21 higher in DRY (1)*
*66 total significant genes*


**These are the BinPacker transcript IDs for the significant genes:**

	detags <- rownames(b)[as.logical(de)]

	detags

 [1] "BINPACKER.10034.2"  "BINPACKER.10141.3"  "BINPACKER.1061.6"   "BINPACKER.10743.2" 
 
 [5] "BINPACKER.11512.1"  "BINPACKER.11560.2"  "BINPACKER.116235.1" "BINPACKER.12709.1"
  
 [9] "BINPACKER.13054.2"  "BINPACKER.13701.1"  "BINPACKER.13806.1"  "BINPACKER.14160.1" 
 
[13] "BINPACKER.147548.1" "BINPACKER.15365.1"  "BINPACKER.15806.1"  "BINPACKER.16191.1" 

[17] "BINPACKER.17022.3"  "BINPACKER.17734.1"  "BINPACKER.17981.2"  "BINPACKER.17992.1"
 
[21] "BINPACKER.1818.1"   "BINPACKER.1846.1"   "BINPACKER.18534.1"  "BINPACKER.18622.1" 

[25] "BINPACKER.20114.1"  "BINPACKER.20530.1"  "BINPACKER.20656.1"  "BINPACKER.20716.2" 

[29] "BINPACKER.21794.1"  "BINPACKER.22521.1"  "BINPACKER.23756.2"  "BINPACKER.23790.1" 

[33] "BINPACKER.24398.1"  "BINPACKER.24914.1"  "BINPACKER.28731.1"  "BINPACKER.28.2" 
   
[37] "BINPACKER.2960.1"   "BINPACKER.31087.1"  "BINPACKER.31815.1"  "BINPACKER.3452.1"  

[41] "BINPACKER.3510.3"   "BINPACKER.35470.1"  "BINPACKER.3957.3"   "BINPACKER.42718.1" 

[45] "BINPACKER.4449.4"   "BINPACKER.4855.1"   "BINPACKER.49203.1"  "BINPACKER.52106.1" 

[49] "BINPACKER.5280.2"   "BINPACKER.56553.1"  "BINPACKER.5662.4"   "BINPACKER.58702.1" 

[53] "BINPACKER.6383.3"   "BINPACKER.6494.2"   "BINPACKER.66588.1"  "BINPACKER.6740.3"  

[57] "BINPACKER.6807.1"   "BINPACKER.724.4"    "BINPACKER.7740.1"   "BINPACKER.87639.1"
 
[61] "BINPACKER.9218.3"   "BINPACKER.93518.1"  "BINPACKER.9604.1"   "BINPACKER.9726.1"  

[65] "BINPACKER.9726.2"   "BINPACKER.9961.2"



**This is a command to print out the 66 significant results:**

	topTags(et, n=66, adjust.method="BH", sort.by="PValue")

*RESULTS:*

Comparison of groups:  2-1 

  logFC ;     logCPM   ;    PValue    ;    FDR
                       
BINPACKER.15365.1  -3.702503  0.04719296 7.687418e-16 5.311007e-11

BINPACKER.2960.1   -4.267653  1.14652643 5.973998e-14 2.063628e-09

BINPACKER.21794.1   2.434172  3.11662976 2.455055e-12 4.413836e-08

BINPACKER.28731.1   2.483909  1.63439651 2.555523e-12 4.413836e-08

BINPACKER.17981.2  -2.975289  0.43560006 4.555071e-12 6.293924e-08

BINPACKER.5662.4    2.061044  2.41894040 1.142498e-11 1.315530e-07

BINPACKER.87639.1   2.682415  0.34504725 1.988020e-11 1.962091e-07

BINPACKER.9961.2   -2.425763  1.99801667 8.681049e-11 7.496845e-07

BINPACKER.3452.1   -2.506611 -0.13960225 4.638913e-10 3.560984e-06

BINPACKER.724.4    -2.162465  2.66706772 1.204519e-09 8.321661e-06

BINPACKER.9604.1   -2.582197  0.54693424 1.253510e-08 7.872838e-05

BINPACKER.31087.1  -2.907820 -0.85824933 1.763985e-08 9.741421e-05

BINPACKER.24398.1  -2.440435 -0.68906346 1.833029e-08 9.741421e-05

BINPACKER.35470.1   2.366595  1.78613511 3.823434e-08 1.886783e-04

BINPACKER.9726.1   -3.473617 -0.10671828 5.169043e-08 2.380758e-04

BINPACKER.9218.3   -1.578145  1.52518777 6.403276e-08 2.764894e-04

BINPACKER.18534.1  -2.332353  1.34565276 1.194036e-07 4.852494e-04

BINPACKER.52106.1   2.096392 -0.54223174 1.778343e-07 6.825577e-04

BINPACKER.17022.3  -2.898614 -0.56116043 2.760203e-07 1.003653e-03

BINPACKER.3957.3    6.309200  1.57909931 2.960271e-07 1.022581e-03

BINPACKER.13806.1  -2.442034 -0.38086864 3.537466e-07 1.132135e-03

BINPACKER.7740.1   -2.790117  1.09464938 3.605159e-07 1.132135e-03

BINPACKER.10034.2  -4.420235  0.38690629 4.109800e-07 1.234495e-03

BINPACKER.11560.2  -1.464633  2.04967724 5.750024e-07 1.655216e-03

BINPACKER.13701.1  -1.312101  1.80418082 8.244510e-07 2.278354e-03

BINPACKER.3510.3   -2.162798  0.90579216 1.108665e-06 2.945937e-03

BINPACKER.15806.1  -1.699571  1.06224055 1.348431e-06 3.393105e-03

BINPACKER.17992.1  -2.542206  0.65319164 1.375178e-06 3.393105e-03

BINPACKER.9726.2   -2.118980  0.56040947 1.462141e-06 3.483273e-03

BINPACKER.116235.1  2.212301  0.30149366 1.712892e-06 3.944619e-03

BINPACKER.6383.3   -2.093451  1.26954525 1.867001e-06 4.160823e-03

BINPACKER.20716.2  -4.203965 -0.56629460 2.665336e-06 5.754376e-03

BINPACKER.20114.1  -1.661498  0.50070312 2.851677e-06 5.970115e-03

BINPACKER.18622.1  -1.645134  1.70433063 3.131289e-06 6.362688e-03

BINPACKER.4449.4    3.428414 -0.53754364 3.416400e-06 6.743682e-03

BINPACKER.24914.1  -2.211028 -0.15945079 5.157351e-06 9.826820e-03

BINPACKER.31815.1  -1.904640 -0.76976547 5.262818e-06 9.826820e-03

BINPACKER.6740.3   -3.089797 -0.43430811 5.738824e-06 1.043363e-02

BINPACKER.28.2      4.183445  2.29531990 5.935987e-06 1.051537e-02

BINPACKER.20530.1  -1.625831  0.54456222 6.488830e-06 1.120734e-02

BINPACKER.20656.1  -1.910151 -0.53088156 7.224770e-06 1.217409e-02

BINPACKER.4855.1   -1.339923  4.02496858 7.570742e-06 1.233782e-02

BINPACKER.1846.1   -3.280068 -0.79214722 7.679106e-06 1.233782e-02

BINPACKER.6494.2   -3.363298  0.02949406 7.995890e-06 1.255482e-02

BINPACKER.56553.1   1.471615  0.17203428 9.517793e-06 1.461235e-02

BINPACKER.93518.1   1.710779 -0.79262563 1.046158e-05 1.571216e-02

BINPACKER.11512.1   1.186807  3.65446885 1.153285e-05 1.695255e-02

BINPACKER.66588.1   1.850643 -0.34747423 1.188372e-05 1.710438e-02

BINPACKER.1818.1   -1.712717  3.28859652 1.439917e-05 2.030195e-02

BINPACKER.10743.2  -1.915392 -0.52483309 1.514198e-05 2.060668e-02

BINPACKER.42718.1   1.541791  0.50677177 1.535060e-05 2.060668e-02

BINPACKER.13054.2  -1.147098  2.69670509 1.551012e-05 2.060668e-02

BINPACKER.6807.1   -1.329986  2.10573616 1.636913e-05 2.133762e-02

BINPACKER.49203.1   1.638561 -0.03463627 1.908001e-05 2.441076e-02

BINPACKER.14160.1  -2.050541  0.60346386 2.275675e-05 2.858537e-02

BINPACKER.147548.1  1.744028 -0.00714171 2.420523e-05 2.986191e-02

BINPACKER.23756.2   1.264987  3.46849716 2.490620e-05 3.011647e-02

BINPACKER.12709.1   3.906196  2.61069936 2.528341e-05 3.011647e-02

BINPACKER.16191.1  -1.431323  0.92602592 2.918057e-05 3.416945e-02

BINPACKER.10141.3  -3.283133 -1.19055421 3.193816e-05 3.677520e-02

BINPACKER.5280.2    3.873777  0.25690923 3.319288e-05 3.759339e-02

BINPACKER.23790.1  -1.755781 -0.27493493 4.043437e-05 4.505628e-02

BINPACKER.22521.1  -1.840899 -0.05571504 4.126304e-05 4.524983e-02

BINPACKER.1061.6   -1.807291  1.94295119 4.631104e-05 4.926979e-02

BINPACKER.58702.1   1.779863 -0.49958790 4.635512e-05 4.926979e-02

BINPACKER.17734.1  -1.659928  2.10923578 4.721537e-05 4.942377e-02


**plot:**

	plotSmear(et, de.tags=detags, main="DE genes, all data, FDR 0.05")
	abline(h = c(-2, 2), col = "blue")
	abline(h = c(-4, 4), col = "blue")

-
***Now that I have the transcript IDs which are significantly differentially expressed between groups, the following will allow me to quickly view side by side the transcript ID, gene ID, & gene info:**

	cat NEWESTmusHITS.txt | awk '{print $1 "\t" $7 "\t" $10 "\t" $11}' > NEWESTGENEmusHITS.txt

	cat NEWESTGENEmusHITS.txt | wc -l

*45636*


	cat NEWESTGENEmusHITS.txt | sort -uk1,1 > NEWESTSORTEDgeneHITS.txt

	cat NEWESTSORTEDgeneHITS.txt | wc -l

*37744*


**Then I can view the gene ID and gene info for each transcript ID that was significant, for example:**

	grep "BINPACKER.9726.2" SORTEDgeneHITS.txt


***I also retrieved the BINPACKER transcripts which did not have a matching Ensembl ID in the DTE:***

**This was done by the following command in the folder that contains my transcriptome assembly, for example:** 

	grep "BINPACKER.3452.1" good.BINPACKER.cdhit.fasta -A 1

>BINPACKER.3452.1
CTCCAAGGGATCGGCTTCAAATGGACCTATAACACAGTGATGCGTCCA....

***All if these grepped sequences and their results are in a different file: DTEno-matchBLASTnSequences.md***

