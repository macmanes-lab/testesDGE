#Manuscript: Differential Gene Expression for TESTES


______

**Part 2) Differential Gene Expression Analysis**

#a. Map All Raw Read Sets from WET and DRY Mice to the Testes Transcriptome



***First I counted all the reads in the following dataset to confirm that the left and right read sets had the same numbers of sequences, and they all did.***

**Example:**

	grep @HWI -c 1357T_R1.fastq

*20603232 sequences*

	grep @HWI -c 1357T_R2.fastq

*20603232 sequences*


***Mapping with SALMON v07.2***

**Index Step only has to be done once for the Testes transcriptome:**


	salmon index -t good.BINPACKER.cdhit.fasta -i NEWTESTES_salmon.idx --type quasi -k 31


**Quantification Step has to be done for each read set:**

	salmon quant -p 32 --seqBias --gcBias -i NEWTESTES_salmon.idx -l a -1 3333T_R1.fastq -2 3333T_R2.fastq -o NEWmaps/salmon/salmon_3333T


	salmon quant -p 32 --seqBias --gcBias -i NEWTESTES_salmon.idx -l a -1 T2322_R1.fastq -2 T2322_R2.fastq -o NEWmaps/salmon/salmon_2322T


	salmon quant -p 32 --seqBias --gcBias -i NEWTESTES_salmon.idx -l a -1 382T_R1.fastq -2 382T_R2.fastq -o NEWmaps/salmon/salmon_382T


	salmon quant -p 32 --seqBias --gcBias -i NEWTESTES_salmon.idx -l a -1 381T_R1.fastq -2 381T_R2.fastq -o NEWmaps/salmon/salmon_381T


	salmon quant -p 32 --seqBias --gcBias -i NEWTESTES_salmon.idx -l a -1 376T_R1.fastq -2 376T_R2.fastq -o NEWmaps/salmon/salmon_376T


	salmon quant -p 32 --seqBias --gcBias -i NEWTESTES_salmon.idx -l a -1 366T_R1.fastq -2 366T_R2.fastq -o NEWmaps/salmon/salmon_366T


	salmon quant -p 32 --seqBias --gcBias -i NEWTESTES_salmon.idx -l a -1 349T_R1.fastq -2 349T_R2.fastq -o NEWmaps/salmon/salmon_349T


	salmon quant -p 32 --seqBias --gcBias -i NEWTESTES_salmon.idx -l a -1 209T_R1.fastq -2 209T_R2.fastq -o NEWmaps/salmon/salmon_209T


	salmon quant -p 32 --seqBias --gcBias -i NEWTESTES_salmon.idx -l a -1 265T_R1.fastq -2 265T_R2.fastq -o NEWmaps/salmon/salmon_265T


	salmon quant -p 32 --seqBias --gcBias -i NEWTESTES_salmon.idx -l a -1 383T_R1.fastq -2 383T_R2.fastq -o NEWmaps/salmon/salmon_383T


	salmon quant -p 32 --seqBias --gcBias -i NEWTESTES_salmon.idx -l a -1 384T_R1.fastq -2 384T_R2.fastq -o NEWmaps/salmon/salmon_384T


	salmon quant -p 32 --seqBias --gcBias -i NEWTESTES_salmon.idx -l a -1 102T_R1.fastq -2 102T_R2.fastq -o NEWmaps/salmon/salmon_102T


	salmon quant -p 32 --seqBias --gcBias -i NEWTESTES_salmon.idx -l a -1 400T_R1.fastq -2 400T_R2.fastq -o NEWmaps/salmon/salmon_400T


	salmon quant -p 32 --seqBias --gcBias -i NEWTESTES_salmon.idx -l a -1 1357T_R1.fastq -2 1357T_R2.fastq -o NEWmaps/salmon/salmon_1357T


	salmon quant -p 32 --seqBias --gcBias -i NEWTESTES_salmon.idx -l a -1 1358T_R1.fastq -2 1358T_R2.fastq -o NEWmaps/salmon/salmon_1358T


	salmon quant -p 32 --seqBias --gcBias -i NEWTESTES_salmon.idx -l a -1 1359T_R1.fastq -2 1359T_R2.fastq -o NEWmaps/salmon/salmon_1359T


	salmon quant -p 32 --seqBias --gcBias -i NEWTESTES_salmon.idx -l a -1 13T_R1.fastq -2 13T_R2.fastq -o NEWmaps/salmon/salmon_13T


	salmon quant -p 32 --seqBias --gcBias -i NEWTESTES_salmon.idx -l a -1 343T_R1.fastq -2 343T_R2.fastq -o NEWmaps/salmon/salmon_343T


	salmon quant -p 32 --seqBias --gcBias -i NEWTESTES_salmon.idx -l a -1 344T_R1.fastq -2 344T_R2.fastq -o NEWmaps/salmon/salmon_344T


	salmon quant -p 32 --seqBias --gcBias -i NEWTESTES_salmon.idx -l a -1 355T_R1.fastq -2 355T_R2.fastq -o NEWmaps/salmon/salmon_355T


	salmon quant -p 32 --seqBias --gcBias -i NEWTESTES_salmon.idx -l a -1 888T_R1.fastq -2 888T_R2.fastq -o NEWmaps/salmon/salmon_888T



#b. Genrate Gene x Transcript ID matrix
	
The next step is to generate a matrix table relating each transcript ID to a gene ID. This is critical for a DGE approach, because multiple transcripts can be transcribed by each gene. Therefore, I have to label each transcript ID with the appropriate gene ID. 

In order to do so, I will download the latest Mus cds file from Ensembl and do a BLASTn of the Mus database to my testes transcriptome. This will assign matching transcript IDs in my testes transcriptome to Transcript IDs and Gene IDs in the Mus transcriptome. I will then select only the Gene IDs from the Mus transcriptome for each of the P. eremicus testes transcript IDs. This will provide me with a Gene x transcript ID table. This table is necessary for tximport.  This table will be used for DGE analysis.


The Ensembl Mus database was dowloaded at 9:15 am EST on 9/5/16
 
Index of /pub/release-85/fasta/mus_musculus/cds
version date: 7/11/16, 3:05:00 AM

Mus_musculus.GRCm38.cds.all.fa.gz
15.9 MB	7/7/16, 10:04:00 AM


	mkdir MUS

	wget ftp://ftp.ensembl.org/pub/release-85/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz

	gzip -d Mus_musculus.GRCm38.cds.all.fa.gz

	makeblastdb -in /mnt/data3/lauren/TESTESdge/MUS/Mus_musculus.GRCm38.cds.all.fa -out MUScds -dbtype nucl

**This is a nucleotide database, and I have a nucleotide query for my transcriptome, so it is a BLASTn search:**

	blastn -query /mnt/data3/lauren/good.BINPACKER.cdhit.fasta \
	-db /mnt/data3/lauren/TESTESdge/MUS/MUScds \
	-max_target_seqs 1 \
	-outfmt '6 qseqid pident evalue stitle' \
	-evalue 1e-5 -num_threads 10 | tee NEWESTmusHITS.txt


**Next I have to modify the BLASTn results to make my matrix table:**

	cat NEWESTmusHITS.txt | wc -l 

*45636*


	cat NEWESTmusHITS.txt | awk '{print $1 "\t" $7}' > NEWESTSHORTmusHITS.txt


	cat NEWESTSHORTmusHITS.txt | wc -l

*45636* 


	cat NEWESTSHORTmusHITS.txt | sort -uk1,1 > NEWESTSORTEDmusHITS.txt


	cat NEWESTSORTEDmusHITS.txt | wc -l

*37744*

	sed 's/gene://g' NEWESTSORTEDmusHITS.txt > NEWESTFinalMUS.txt


	cat NEWESTFinalMUS.txt | wc -l

*37744*


***NEWESTFinalMUS.txt is file that only has one transcriptID per MUS gene ID that has matched; therefore this is my Matrix file.***

*There are 37,744 matches between transcript IDs in my testes transcriptome and the matches to the MUS transcriptome.*


**Calculation for number of Mus gene matches:**

	cat NEWESTFinalMUS.txt | awk '{print $2}'| sort -u | wc -l 

*14218*

Of these matches 37,744 matches to transcript IDs in the testes transcriptome, a little less than a third of them (14,218) are for unique gene ID matches to MUS ( representing ~14k gene IDs) 

