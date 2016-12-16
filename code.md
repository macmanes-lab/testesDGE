#Manuscript: Differential Gene Expression for TESTES

**Part 1) Assembly of a Testes Transcriptome**

**Part 2) Differential Gene Expression Analysis**


-
**location of raw testes reads:**

lauren@davinci:/mnt/data3/lauren/TESTESdge/testes.R1.fastq

lauren@davinci:/mnt/data3/lauren/TESTESdge/testes.R2.fastq

*READ SET INFORMATION*

There are 45759144 total fragments in R1 AND there are 45759144 total fragments in R2

So it is a ~45.7 million paired-end read data set

Mapping rate (SALMON) of raw reads to final un-annotated testes transcriptome: 85.4604%


**location of final un-annotated testes transcriptome:**

lauren@davinci:/mnt/data3/lauren/TESTESassemblyNOV16/good.BINPACKER.cdhit.fasta 

*ASSEMBLY INFORMATION*

Number of Transcripts: 155134

Transrate Score: 0.3352 

BUSCO metrics: 77% SCO ;	27 % DCO ;	5.9% Fragmented ; 16% Missing


**location of final dammit-annotated testes transcriptome:**

lauren@Pinky:/home/lauren/Documents/Nov2016Dammit/dammit/good.BINPACKER.cdhit.fasta.dammit/good.BINPACKER.cdhit.fasta.dammit.fasta

*ASSEMBLY INFORMATION*

good.BINPACKER.cdhit.fasta.dammit.fasta 
Number of transcripts: 155134


-
***Part 1) Assembly of a Testes Transcriptome***


**The testes transcriptome was assembled several ways to determine an optimal assembly**


***Each time a new assembly was built or an existing assembly was optimized, we generated statistics and parameters to evaluate the relative quality and completeness of the transcriptome.  I have included code once for this first assembly as an example for these indices which we generated for all assembly versions.***

	a. Reads	

Error Correction in Rcorrector

Rcorrector (v1.0.1)
perl /share/Rcorrector/run_rcorrector.pl -k 31 -t 30 \
-1 testes.R1.fastq \
-2 testes.R2.fastq
	
Error Corrected Reads: 
testes.R1.cor.fq
testes.R2.cor.fq


	b. Trinity Assembly

Trinity (v.2.2.0)

Trinity --SS_lib_type RF --seqType fq --max_memory 40G --trimmomatic --CPU 30 --full_cleanup --output Testes_trinity \
--left testes.R1.cor.fq \
--right testes.R2.cor.fq \
--quality_trimming_params "ILLUMINACLIP:/share/trinityrnaseq/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25"


Name of Trinity Original Assembly: Testes_trinity.Trinity.fasta

abyss-fac Testes_trinity.Trinity.fasta

n		n:500	L50		min		N80		N50		N20  E-size max   sum     name
397557	127542	26874	500		868		1963	4036 2658   19259 185.6e6  Testes_trinity.Trinity.fasta


number of transcripts: 397557	

BUSCO (v1.1b1) for original assmbly:
python3 /share/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 16 -l /share/BUSCO_v1.1b1/vertebrata \
-o TestesTrinity_busco -in /mnt/data3/lauren/NEWtranscriptome/Testes_trinity.Trinity.fasta

Summarized benchmarks in BUSCO notation:
	C:79%[D:44%],F:5.3%,M:14%,n:3023
	
TRANSRATE (transrate v1.0.1)
/share/transrate-1.0.1-linux-x86_64/transrate --output NEWEST_TESTES_transrate -t 16 \
--assembly Testes_trinity.Trinity.fasta \
--left testes.R1.cor.fq \
--right testes.R2.cor.fq 

TRANSRATE ASSEMBLY SCORE     0.196

TRANSRATE OPTIMAL SCORE      0.3329


Name of Transrate Optimized Assembly: good.Testes_trinity.Trinity.fasta

abyss-fac NEWEST_TESTES_transrate/Testes_trinity.Trinity/good.Testes_trinity.Trinity.fasta

n		n:500	L50		min		N80		N50		N20  E-size max   sum     name
342153	114100	24774	500		831		1828	3731 2457   19259 158.7e6 good.Testes_trinity.Trinity.fasta

number of transcripts: 342153

BUSCO on transrate corrected transcriptome:
python3 /share/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 16 -l /share/BUSCO_v1.1b1/vertebrata \
-o NewTestesTransrate_busco -in /mnt/data3/lauren/NEWtranscriptome/NEWEST_TESTES_transrate/Testes_trinity.Trinity/good.Testes_trinity.Trinity.fasta

Summarized benchmarks in BUSCO notation:
        C:79%[D:43%],F:5.3%,M:14%,n:3023


	c. BinPacker Assembly

SKEWER: skewer (0.1.127) 

skewer -l 25 -m pe -o skewer --mean-quality 2 --end-quality 2 -t 30 \
-x /share/trinityrnaseq/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE-2.fa \
testes.R1.cor.fq testes.R2.cor.fq


BINPACKER: BinPacker (v1.0)

/opt/BinPacker/BinPacker -d -q -s fq -p pair -m RF -k 25 -g 200 -o TESETES_binpacker \
-l skewer-trimmed-pair1.fastq \
-r skewer-trimmed-pair2.fastq

Name of Assembly: BinPacker.fa

number of transcripts: 170958	


	d. Transfuse Assembly

TRANSFUSE (v.0.5.0):

/opt/transfuse-0.5.0-linux-x86_64/transfuse -t 16 -i 0.98 -o transfuse_TESTES \
-l skewer-trimmed-pair1.fastq \
-r skewer-trimmed-pair2.fastq \
-a /home/lauren/Documents/TESTESannotation/TESETES_binpacker/BinPacker.fa,/home/lauren/Documents/TESTESannotation/Testes_trinity.Trinity.fasta


Name of Assembly: transfuse_TESTES_cons.fa

Number of Transcripts: 361796


	e. Additional optimization methods
	


***CD-hit-est***

Example for CD-hit-est on the original binpacker transrate corrected assembly: goodBINPACKER.fa

cd-hit-est -M 5000 -T 14 -c .97 -i goodBINPACKER.fa -o BINPACKER.cdhit.fasta

Assembly name: BINPACKER.cdhit.fasta

Number of Transcripts: 155588	


***Transrate correction***

/opt/transrate-1.0.4beta/transrate -o transrate_binpackerCDhit -t 32 \
-a /home/lauren/Documents/TESTESannotation/FIXING/BINPACKER.cdhit.fasta \
--left /home/lauren/Documents/TESTESannotation/FIXING/skewer-trimmed-pair1.fastq --right /home/lauren/Documents/TESTESannotation/FIXING/skewer-trimmed-pair2.fastq

The most important stats for the transrate run are below: 
transrate assembly score (0.3317)
optimal score (0.3352)

abyss-fac good.BINPACKER.cdhit.fasta

n       n:500   L50     min     N80     N50     N20     E-size  max     sum     
155134  94099   20181   500     876     1927    3953    2595    20787   136.7e6


***BUSCO***

python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 32 -l /opt/BUSCO_v1.1b1/vertebrata \
-o BUSCO_TRANSRATEbinpackerCDhit -in /home/lauren/Documents/TESTESannotation/FIXING/transrate_binpackerCDhit/BINPACKER.cdhit/good.BINPACKER.cdhit.fasta


BUSCO results:
Summarized benchmarks in BUSCO notation:
        C:77%[D:27%],F:5.9%,M:16%,n:3023

Representing:
        1524    Complete Single-copy BUSCOs
        819     Complete Duplicated BUSCOs
        180     Fragmented BUSCOs
        500     Missing BUSCOs
        3023    Total BUSCO groups searched



This is where the new transcriptome assembly was copied from in PINKY:

lauren@davinci:/mnt/data3/lauren/TESTESassemblyNOV16/FIXING/transrate_binpackerCDhit/BINPACKER.cdhit/
good.BINPACKER.cdhit.fasta

It is now here on DAVINCI:
lauren@davinci:/mnt/data3/lauren/TESTESassemblyNOV16/good.BINPACKER.cdhit.fasta

And here on DAVINCI:
lauren@davinci:/mnt/data3/lauren/good.BINPACKER.cdhit.fasta



***Filtering by Expression Levels***

Example of filtering by TPM values of 0.5 or 1.0 for the trinity assembly

KALLISTO (v0.42.4):

kallisto index -i Testes.idx Testes_trinity.Trinity.fasta

kallisto quant -t 32 -i Testes.idx -o TestesOnly.output testes.R1.cor.fq testes.R2.cor.fq 


SALMON (v0.6.0):

/share/salmon-0.6.0/build/src/salmon index -t Testes_trinity.Trinity.fasta -i Testes_salmon.idx --type quasi -k 31

/share/salmon-0.6.0/build/src/salmon quant -p 32 -i Testes_salmon.idx -l IU -1 testes.R1.cor.fq -2 testes.R2.cor.fq -o salmon_TestesOnly


*Filtering at TPM 0.5*


sh FILE (The output of the script is TestesOnlyHALF.trinity.Trinity.fasta)
Filtering to keep only transcripts with TPM>0.5 		
file name: NEWTESTEStpmHALF.sh
sh ./NEWTESTEStpmHALF.sh 

-

awk '0.5>$5{next}1' TestesOnly.output/abundance.tsv | awk '{print $1}' > NEWkallist
awk '0.5>$4{next}1' salmon_TestesOnly/quant.sf | sed  '1,10d' | awk '{print $1}' > NEWsalist
cat NEWkallist NEWsalist | sort -u > NEWuniq_list
sed -i ':begin;N;/[ACTGNn-]\n[ACTGNn-]/s/\n//;tbegin;P;D' Testes_trinity.Trinity.fasta

for i in $(cat NEWuniq_list);
   do grep --no-group-separator --max-count=1 -A1 -w $i Testes_trinity.Trinity.fasta >> NEWTestesOnlyHALF.trinity.Trinity.fasta;
done

-

Name of Assembly: NEWTestesOnlyHALF.trinity.Trinity.fasta

Number of Transcripts: 280298


*Filtering at TPM 1.0*


sh FILE (The output of the script is TestesOnlyONE.trinity.Trinity.fasta)
Filtering to keep only transcripts with TPM>1 		
file name: NEWTESTEStpmONE.sh
sh ./NEWTESTEStpmONE.sh 

-

awk '1>$5{next}1' TestesOnly.output/abundance.tsv | awk '{print $1}' > kallistONE
awk '1>$4{next}1' salmon_TestesOnly/quant.sf | sed  '1,10d' | awk '{print $1}' > salistONE
cat kallistONE salistONE | sort -u > uniq_listONE
sed -i ':begin;N;/[ACTGNn-]\n[ACTGNn-]/s/\n//;tbegin;P;D' Testes_trinity.Trinity.fasta

for i in $(cat uniq_listONE);
   do grep --no-group-separator --max-count=1 -A1 -w $i Testes_trinity.Trinity.fasta >> NEWTestesOnlyONE.trinity.Trinity.fasta;
done

-

Name of Assembly: NEWTestesOnlyONE.trinity.Trinity.fasta

Number of Transcripts: 122825	




***The chosen TESTES assembly was generated and optimized in the following ways:***

1) Raw reads were error corrected with R-corrector

2) Error corrected reads were trimmed with skewer

3) Assembly was done with BinPacker

4) BinPacker Assembly was optimized with Transrate to retain highly supported contigs

5) Transrate Corrected assembly was put through CD-hit-est to reduce redundancies

6) This assembly was optimized with Transrate to retain highly supported contigs


This is where the new transcriptome assembly was copied from in Davinci:
lauren@davinci:/mnt/data3/lauren/TESTESassemblyNOV16/FIXING/transrate_binpackerCDhit/BINPACKER.cdhit/
good.BINPACKER.cdhit.fasta

It is now in two places in DAVINCI:
lauren@davinci:/mnt/data3/lauren/TESTESassemblyNOV16/good.BINPACKER.cdhit.fasta
lauren@davinci:/mnt/data3/lauren/good.BINPACKER.cdhit.fasta

***Assembly Name: good.BINPACKER.cdhit.fasta***

Number of Transcripts: 155134

Transrate Score: 0.3352 

BUSCO metrics: 77% SCO ;	27 % DCO ;	5.9% Fragmented ; 16% Missing


	f. Mapping Rates for Chosen Optimal Assembly

***The reads were mapped back to this chosen assembly with Salmon to determine the mapping rate as a final quality check:***

Salmon v0.5.1 mapping of reads to assembly

/opt/salmon-0.5.1/bin/salmon index -t /home/lauren/Documents/TESTESannotation/FIXING/transrate_binpackerCDhit/BINPACKER.cdhit/good.BINPACKER.cdhit.fasta -i goodBINPACKERcdhit_salmon.idx --type quasi -k 31

/opt/salmon-0.5.1/bin/salmon quant -p 32 -i goodBINPACKERcdhit_salmon.idx -l IU -1 skewer-trimmed-pair1.fastq -2 skewer-trimmed-pair2.fastq -o salmon_goodBINPACKERcdhit

The % mapping results of salmon are high: 

[2016-11-08 19:52:43.335] [jointLog] [info] Mapping rate = 92.1438%



	g. ASSEMBLY ANNOTATION with DAMMIT
	
*The transcriptome was moved to Pinky for the Dammit assembly:*

scp lauren@davinci:/mnt/data3/lauren/TESTESassemblyNOV16/good.BINPACKER.cdhit.fasta .


*Therefore the transcriptome is here on Pinky:*

lauren@Pinky:~/Documents/Nov2016Dammit/good.BINPACKER.cdhit.fasta


**The dammit directory is here:**

mkdir dammit
cd dammit

**Install the dammit databases:**

dammit databases --install --database-dir /home/lauren/Documents/Nov2016Dammit/dammit --full --busco-group vertebrata 


**Annotate transcriptome:**

dammit annotate /home/lauren/Documents/Nov2016Dammit/good.BINPACKER.cdhit.fasta --busco-group vertebrata --n_threads 36 --database-dir /home/lauren/Documents/Nov2016Dammit/dammit/ --full


*output is new annotated assembly check destination:*

/home/lauren/Documents/Nov2016Dammit/dammit/good.BINPACKER.cdhit.fasta.dammit/good.BINPACKER.cdhit.fasta.dammit.fasta


abyss-fac good.BINPACKER.cdhit.fasta.dammit.fasta 

n	n:500	L50	min	N80	N50	N20	E-size	max	sum	name
155134	94099	20181	500	876	1927	3953	2595	20787	136.7e6	good.BINPACKER.cdhit.fasta.dammit.fasta


***I counted the ORFs meeting the minimum critical length in the pep file:***
  
grep ">" good.BINPACKER.cdhit.fasta.transdecoder.pep | wc -l

75482 transcripts have an orf (48.7% is 75482/155134)


grep "complete" good.BINPACKER.cdhit.fasta.transdecoder.pep | wc -l


43028 transcripts have a complete orf (57.0% is 43028/75482)


***I counted the pfam results:***

sort -uk1,1 good.BINPACKER.cdhit.fasta.pfam.csv.gff3 | wc -l

25675 transcripts have a pfam hit (16.6% is 25675/155134)


***I counted the LAST (akin to blast) results:***

sort -uk1,1 good.BINPACKER.cdhit.fasta.x.uniref.maf.best.csv.gff3 | wc -l

62865 transcripts have uniref90 hits (40.5% is 62865/155134)


***I counted the rfam results:***
sort -uk1,1 good.BINPACKER.cdhit.fasta.rfam.tbl.gff3 | wc -l

937 transcripts have rfam hits (ncRNAs) (0.6% is 937/155134)


***I counted the OrthoDB results:***

sort -uk1,1 good.BINPACKER.cdhit.fasta.x.orthodb.maf.gff3 | wc -l

51806

51806 transcripts have ortodb hits (33.4% is 51806/155134)

***I also recorded the BUSCO results (but they are the same as the un-annotated :***

testes transcriptome BUSCO results)


***I counted the overall annotated results for dammit:*** 

sort -uk1,1 good.BINPACKER.cdhit.fasta.dammit.gff3 | wc -l

77915 transcripts have dammit annotation hits of some kind from the five DBs (50.2% is 77915/155134)


______

**Part 2) Differential Gene Expression Analysis**

	a. Map All Raw Read Sets from WET and DRY Mice to the Testes Transcriptome

*Note that the testes read set used to build the assembly was from mouse number 335 (wet desert conditions)*

***First I counted all the reads in the following dataset to confirm that the left and right read sets had the same numbers of sequences, and they all did.***

example command:

grep @HWI -c 1357T_R1.fastq

20603232

grep @HWI -c 1357T_R2.fastq

20603232


***Raw Read Data Sets***
	
3333T_R1.fastq

3333T_R2.fastq
	

T2322_R2.fastq

T2322_R1.fastq


382T_R2.fastq

382T_R1.fastq
	

381T_R2.fastq

381T_R1.fastq


376T_R2.fastq

376T_R1.fastq
	

366T_R2.fastq

366T_R1.fastq


349T_R2.fastq

349T_R1.fastq


209T_R1.fastq

209T_R2.fastq
	

265T_R1.fastq

265T_R2.fastq


383T_R1.fastq

383T_R2.fastq


384T_R1.fastq

384T_R2.fastq


102T_R1.fastq

102T_R2.fastq

	
400T_R1.fastq

400T_R2.fastq


1357T_R1.fastq

1357T_R2.fastq


1358T_R1.fastq

1358T_R2.fastq


1359T_R1.fastq

1359T_R2.fastq


13T_R1.fastq

13T_R2.fastq


343T_R1.fastq

343T_R2.fastq


344T_R1.fastq

344T_R2.fastq


355T_R1.fastq

355T_R2.fastq


888T_R1.fastq

888T_R2.fastq


999T_R1.fastq

999T_R2.fastq



***Mapping with SALMON (v0.6.0) THEN replaced with new salmon v07.2:***

*Index Step only has to be done once for the Testes transcriptome:*


I also mapped the original read set with salmon, for a quality check comparison:


salmon index -t good.BINPACKER.cdhit.fasta -i NEWTESTES_salmon.idx --type quasi -k 31

salmon quant -p 32 --seqBias --gcBias -i NEWTESTES_salmon.idx -l a -1 testes.R1.fastq -2 testes.R2.fastq -o NEWmaps/salmon_testesR1R2

*the output files are located in this location:*

lauren@davinci:/mnt/data3/lauren/TESTESdge/NEWmaps/salmon


*Quantification Step has to be done for each read set:*

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




	b. Genrate Gene x Transcript ID table
	
The next step is to generate a table relating each transcript ID to a gene ID. This is critical for a DGE approach, because multiple transcripts can be transcribed by each gene. Therefore, I have to label each transcript ID with the appropriate gene ID. 

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

This is a nucleotide database, and I have a nucleotide query for my transcriptome, so it is a blastn


blastn -query /mnt/data3/lauren/good.BINPACKER.cdhit.fasta \
-db /mnt/data3/lauren/TESTESdge/MUS/MUScds \
-max_target_seqs 1 \
-outfmt '6 qseqid pident evalue stitle' \
-evalue 1e-5 -num_threads 10 | tee NEWESTmusHITS.txt



cat NEWESTmusHITS.txt | wc -l 

* 45636


cat NEWESTmusHITS.txt | awk '{print $1}'| wc -l 

* 45636


cat NEWESTmusHITS.txt | awk '{print $1}'| sort -u | wc -l 

* 37744


cat NEWESTmusHITS.txt | awk '{print $1 "\t" $7}' > NEWESTSHORTmusHITS.txt


cat NEWESTSHORTmusHITS.txt | wc -l

* 45636 


cat NEWESTSHORTmusHITS.txt | sort -uk1,1 > NEWESTSORTEDmusHITS.txt


cat NEWESTSORTEDmusHITS.txt | wc -l

* 37744

sed 's/gene://g' NEWESTSORTEDmusHITS.txt > NEWESTFinalMUS.txt

cat NEWESTFinalMUS.txt | wc -l

* 37744


NEWESTFinalMUS.txt is file that only has one transcriptID per MUS gene ID that has matched.

There are 37,744 matches between transcript IDs in my testes transcriptome and the matches to the MUS transcriptome, which is a lot less than the 155,134 total unique transcript IDs in the testes transcriptome, which means many of the testes transcript IDs had no match to the MUS transcriptome.



cat NEWESTFinalMUS.txt | awk '{print $2}'| sort -u | wc -l 

* 14218

Of these matches 37,744 matches to transcript IDs in the testes transcriptome, a little less than a third of them (14,218) are for unique gene ID matches to MUS (therefore, representing ~14k gene IDs) 

******* 

	c. Differential Gene Expression (DGE) with tximport and edgeR:


***DEPENDENCIES:***

source("https://bioconductor.org/biocLite.R")

biocLite("tximport")

library(tximport)

install.packages("readr")

library(readr)

biocLite("edgeR")

library(edgeR)


***import the gene ID matrix***

tx2gene <- read.delim("~/Desktop/NEWSALMONdge/NEWESTFinalMUS.csv", header=FALSE)
head(tx2gene)


***SET DIRECTORY***

dir <- setwd("/Users/lauren/Desktop/NEWSALMONdge/test")
getwd()
file.path(dir, "salmon", "quant.sf")
dir <- setwd("/Users/lauren/Desktop/NEWSALMONdge/test")
getwd()
file.path(dir, "salmon", "quant.sf")

* This is what should be returned by these commands:
      "/Users/lauren/Desktop/NEWSALMONdge/test/salmon/quant.sf"  
      

***SET samples in downloaded order***


samp <- c("salmon_13T", "salmon_102T", "salmon_209T", "salmon_265T", "salmon_343T", "salmon_344T", "salmon_349T", "salmon_355T", "salmon_366T", "salmon_376T", "salmon_381T", "salmon_382T", "salmon_383T", "salmon_384T", "salmon_400T", "salmon_888T","salmon_999T", "salmon_1357T", "salmon_1358T", "salmon_1359T", "salmon_2322T", "salmon_3333T") 
samp


file.path(dir, "salmon", samp, "quant.sf")
files <- file.path(dir, "salmon", samp, "quant.sf")
files

names(files) <- samp
names(files)


***NOW read in files with tximport***

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, read = read_tsv)
head(txi.salmon$counts)


***Quantification with edgeR***

cts <- txi.salmon$counts
normMat <- txi.salmon$length
normMat <- normMat/exp(rowMeans(log(normMat)))


o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
y <- DGEList(cts)
y$offset <- t(t(log(normMat)) + o)


***assign treatment groups***



group <- c(2,1,2,2,2,1,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1)
group

y <- DGEList(counts=y, group=group)
y

***Scaling plot***

plotMDS.DGEList(y , main = "MDS Plot for Count Data", labels = colnames(y))


* VIEW RESULTS 


***Estimate Common Dispersion***

cds <- estimateCommonDisp(y)
names(cds)

cds$common.dispersion

* 0.5330103  (this is the value returned for common.dispersion)




***Estimate Tagwise Dispersions***
* n integer value should be the nearest integer to this eqn: 50/(#samples - #groups) = 50/(22-2) = 50/20 ~ 2

cds <- estimateTagwiseDisp(cds, prior.n=2)
names(cds)
summary(cds$tagwise.dispersion)

* Values returned: 

  Min. 	1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.1362  0.4216  0.5196  0.5302  0.6228  2.8920 




***I tested whether doing the common & tagwise dispersions a different way effected the results, and it did not:***


b <- estimateCommonDisp(y)
b <- estimateTagwiseDisp(b)
et <- exactTest(b)
topTags(et)

detags <- rownames(topTags(et, n=20))
cpm(b)[detags,]
summary(de <- decideTestsDGE(et, p=0.05, adjust="BH"))
detags <- rownames(b)[as.logical(de)]



***Mean Variance Plot***

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

***Exact Test simplified version***
et <- exactTest(cds, pair=c("1","2"))


summary(de <- decideTestsDGE(et, p=0.05, adjust="BH"))
	 

* RESULTS:

   [,1] 
-1     7
0  14203
1      8

* There are 15 Statistically significant genes
* this means that 7 are expressed more highly in wet	 
* this means that 8 are more highly expressed in dry



* Because there are 15 significant results, just ask for all of them

topTags(et, n=15)


* Results: 

When my fold change is positive it is more highly expressed in dry (dry=2, wet=1)

topTags(et, n=15)

Comparison of groups:  2-1 

                      logFC     logCPM       PValue          FDR
                      
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


* I looked up the ensembl IDs for these significant results, and they are below:

ENSMUSG00000079019.2  Insl3   Insulin-like 3
	
ENSMUSG00000001768.15   Rin2    Ras and Rab Interactor 2
	
ENSMUSG00000054200.6  Ffar4   Free Fatty Acid Receptor 4
	
ENSMUSG00000025479.9    Cyp2e1  cytochrome P450, family 2, subfamily e, polypeptide 1
	
ENSMUSG00000026435.15  Slc45a3 solute carrier family 45, member 3

ENSMUSG00000025020.11 	Slit1   Slit homolog 1 (drosophila)
	
ENSMUSG00000020427.11   Igfbp3  Insulin-like growth factor binding protein 3
	
ENSMUSG00000019997.11    Ctgf    Connective Tissue growth factor
	
ENSMUSG00000031170.14  Slc38a5 solute carrier family 38, member 5
	
ENSMUSG00000040170.13   Fmo2    flavin containing monooxygenase 2 
	
ENSMUSG00000023915.4    Tnfrsf21 tumor necrosis factor receptor superfamily, member 21
	
ENSMUSG00000052974.8    Cyp2f2  cytochrome P450, family 2, subfamily f, polypeptide 2 
	
ENSMUSG00000030830.18   Itgal   integrin alpha L
	
ENSMUSG00000027901.12   Dennd2d DENN/MADD domain containing 2D
	
ENSMUSG00000032554.15   Trf     Transferrin


***plots***

detags <- rownames(cds)[as.logical(de)]
plotSmear(cds, de.tags=detags, ylim=c(-6,6), xlim=c(0,20), frame.plot="false")
abline(h=2, col='blue')
abline(h=-2, col='blue')



__________

**Analysis by BoxPlot of CPM values for differences by TRT in nine genes of interest:**



ENSMUSG00000079019.2  Insl3   Insulin-like 3
	
ENSMUSG00000001768.15   Rin2    Ras and Rab Interactor 2
	
ENSMUSG00000054200.6  Ffar4   Free Fatty Acid Receptor 4
	
ENSMUSG00000026435.15  Slc45a3 solute carrier family 45, member 3
	
ENSMUSG00000020427.11   Igfbp3  Insulin-like growth factor binding protein 3
	
ENSMUSG00000019997.11    Ctgf    Connective Tissue growth factor
	
ENSMUSG00000031170.14  Slc38a5 solute carrier family 38, member 5
	
ENSMUSG00000030830.18   Itgal   integrin alpha L
	
ENSMUSG00000032554.15   Trf     Transferrin




***I graphed all nine of these genes (their cpm values for WET vs. DRY):***

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



*All of these graphs have been exported as PDFs to the DGEmanuscriptCurrent folder (CpmGraphs)*

________

**I also evaluated the relative CPM values of WET vs. DRY for luteinizing hormone**

***But first, I had to find the transcript corresponding to this gene, which I did via BLAST:***

*First I copied my transcriptome (un-annotated testes) into a new folder to do a blast search with it:*

I moved into this folder:

lauren@davinci:/mnt/data3/lauren/LHsearch/db
 
cp /mnt/data3/lauren/TESTESassemblyNOV16/good.BINPACKER.cdhit.fasta .

***I made a BLAST database out of my transcriptome:***

makeblastdb -in /mnt/data3/lauren/LHsearch/db/good.BINPACKER.cdhit.fasta -out TESTEScds -dbtype nucl


This is a nucleotide database, and I have nucleotide queries for my transcriptome, so it is a BLASTn search

I downloaded the Lhgcr and Lhb sequences from NCBI (MUS) to search for these seqs in my Transcriptome


**BLASTn Search for Lhcgr.fasta:**

blastn -query /mnt/data3/lauren/LHsearch/Lhcgr.fasta \
-db /mnt/data3/lauren/LHsearch/db/TESTEScds \
-max_target_seqs 1 \
-outfmt '6 qseqid pident evalue stitle' \
-evalue 1e-5 -num_threads 10 | tee testesLHhits.txt

*RESULT:*
gi|372099093:c88792059-88716861 76.94   0.0     BINPACKER.19.1



**BLASTn Search for Lhb.fasta:**

blastn -query /mnt/data3/lauren/LHsearch/Lhb.fasta \
-db /mnt/data3/lauren/LHsearch/db/TESTEScds \
-max_target_seqs 1 \
-outfmt '6 qseqid pident evalue stitle' \
-evalue 1e-5 -num_threads 10 | tee testesLHhits.txt

*RESULT:*
gi|372099103:45417473-45421855  84.64   1e-177  BINPACKER.9392.2
gi|372099103:45417473-45421855  92.89   7e-100  BINPACKER.9392.2



**I pulled out the sequences from my Testes transcriptome that came up as matches for these two genes:**


***Lhgcr search***

grep "BINPACKER.19.1" good.BINPACKER.cdhit.fasta -A 1

>BINPACKER.19.1
AGATGATCCTGAGGTCAACAGCGGGAAGCCATCCCCGGACATATTCCTCTCCTGTGCACGGAGGTTCAGTCCTGCTCCTGCCCCAGACATGTGCCTCGTCTTTGAAGATGCTCCCAATGGAGTGGAGGCAGCTCTGGCAGCTGGGATGCAGGTTGTCATGGTTCCCGATGAAAACCTGAGTCGAGACTTAACAAGAAAGGCCACAGTGGTGCTGAGTTCCCTGGAGGACTTCCAGCCTGAGCTGTTTGGTCTGCCTGCCTATGAGTGATGGGTGGCAGCCTTGGTGTTGCCATGGGCCTTGTGGCATCCTGAGGGGGGGGCACATGGCATACTGGAGGAGCCACATGGCACCCTTGAGGGGCCACATAGCTCCTGGGGGACACATGGCATCCTGGGGGGCCATATTCAGGGTCCATGCCATTAAGGAGAAAGGGGAAGGCAATGCATAGCACCAGACTGAGCCTGCTTATACACTGTAGCCTGTGAACATGGCCATCCCCGTCTGTGTGGGCTTGGCCATCTCTCTCTGGATTGTGGTCGGGTCAGTAGACAGCAAACAGAGATGGCCCAGAACACAGCTTGAGTACAGTCATGTAATAAAGTGATGTGGTCTGTGAAGAAAATGAACACACACACACACACACACACACACACACACACACACACACAAATTATAATTGCTTAAAAGCCCAGTTCCTATGAAAATGATCAATTGTTCCCCTTTGAAAGATCCTGTCAGGACTTCGGTTGCCCGGCCACAGCTACTGAGCAAGTGTTGTGTGTGAACAAGGAAGGTCATGGTGATTGGAAACAGATTTACCTCCTGAGATTTTTGTCATTCAGTGATACTTTCCTGTCACTTCATATAGGGTCTTCTCCTTTTTTAATAGCAGGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGCGAGCGAGAGAGAGAGAGAGAGAGGAGAGAGAGAGAGAGAGAGGAGAGCACCCATCCTCCCCCGGACATCCCTTCTGAACAAGAGTTCCCATCCTGCCCCGGACTTCCAATTTGGACAAGAGAGCTCCCATCTGGACAAGAGAGAGAGACTTCCTGAATCTGTCAGCTCTGTCTGAACCAAGTGTGTGGATAAGGCCAAGAACGAACCACAAGGAGATGGGCAGACGTCAAGGCAGAAACACATACAACAAAATGAATAGTAATACACCATCACCAGGCCCTAGCCCTCCTCCAACACCTAGACCTGAACATCAGAAATTGGAAGAAGCAGAAGAAAATAGCCTTATGAATGTCATCATGAAGAAGCTAGAGGCTCGTGTAGAGGAAAAGACAAAAAAATGTGAAGAACGCTGTAAACAACTAGAGGAAAGGGCAAACAAATTAGAAGAAATCAAAAAAGTCCTGGAAGAGAACAATAAAATACTGAAAGAAAATCAAGAAAAATCAATGAAACAAATGAAGGAAACAGTCCAAGACCTGAAAAGGGAAATAGAAAAAATGAAGAAGACACAAACAGAGGGAATGCTGGAAATAGAAAATCTGAGAAAACGATTGGGAACTTCAGATGCAAGTATAATCAACAGAATGCAAGAGATGGAAGAGAGAATCTCTAGCGTTGAAGATACAATAGAAGAAATAGATTCATCAGTCAAAGAAAACACTAAAGTCAACAAAGTCATGAACCAAAATGTCCAAGAAATTTGGGACACCATGAAAAGACCAAACCTACGAATAATAGGGGTAGAAGAAGGAGAAGAATACCAACTCAAGGGCACAGAAAATATATTCAACAAGATCATAGAAGAAAACTTTCCCAACTTAAAGAAGGAAATGCCTATGAAGATACAGGAAGCCTATAGAACACCAAACAGACTAGACCCCCAAAAAAAGTCCCCTCGCCACATAATAATTAAACAATTAAACGTACAGAATAAAGAAAGAATATTAAGAGCAGCAAAGGAAAAAGGCCAAGTGACATATAAAGGCAAACCTATCAGAATAACACCCGATTTCTCAATGGAGACTTTGAAAGCCAGAAGGTCCTGGACAGATGTAATGCAGACACTAAGAGACCATGGATGTCAGCCTAGACTAATATACCCAGCAAAACTTTCAATCATCATAGATGGAGTGAACAAGACATTCCATGACAAAGCCAGATTTAAACAATATTTATCCACAAATCCAGCCCTACAGAAAGCACTAGAAGGAAAATTCCAACCTAAGGAAGTCAGATACACCCTCGAAAACACAGGCAATAGATAAAGCCACAGCAGTAAACCCCAACGAAGAAAAGTACACACACATCACCACCAAAAAATAACAGGAATGAACAATCACTGGACATTAATATCCCTCAATATCAATGGACTTAATTCACCTATAAAAAGACATAGGCTTACAGAATGGATACGAAAGCAGGACCCATCTTTCTGCTGCATACAAGAAACACATCTCAAATTCAAAGATAGACACTACCTAAAAATAAAAGGCTGGGAAAAGACTTTCCAATCAAACGGTCTTAAGAAACAAGCGGGTGTAGCCATCCTGATATCCAGCAAAATAGACTTCAAACTAAAATCAATCAAAAGAGATCAAGAAGGGCATTACATACTCATCACAGGAAAGATCCACCAAGATGAAGTCTCAATTCTGAACATTTATGCCCCAAACACAAGGGCACCCACATATGTAAAAGAAACATTATTAAAGCTTAAATCACATATAAAACCCCACACATTAATAGTGGGAGACCTCAACACCCCACTTTCACCACTGGACAGATCCCCCAAATCGAAACTTAACAGAGAAATAAAGGACTTAACTGATGTCATGACTCAAATGGACTTAATCGACATCTACAGAACATTCCATCCTAACAAAAAAGAATATACCTTCTTCTCAGCACCCCATGGAACCTTCTCTAAAATCGACCACATACTTGGTCACAAAACAAATCTAAACAGATACAAAACAATTGGAATAACCTCCTGTGTTCTATCAGACCACCATGGTCTAAAGTTGCATTTCAACAACAACAAAAACTACAGAAAACCTACAATCTCATGGAAACTGAACAATACCCAACTGAATCACCAATGGGTTAAGGAAGAAATAAAGAAAGAAATTAAAGACTTCCTAGAGATCAACGAAAATGAAGACACCACATATCCAAACCTATGGGACACTATGAAAGCAGTACTAAGAGGGAAATTCATAGCACTAAACGCCCACATAAATAAGCTGGAGAAATCTCACACTAGTGACTTAACAGCACACCTGAAAGTTCTAGAACAGGAAGAAGCAAAGTCTCCCAGGAAAAATAGATGCCAGGAAATTATCAAAGTGAGAGCTGAAATCAATAAAATAGAAACAAAGAGAACAATACAAAAAATTAATGAAACAAAGAGTTGGTTCTTTGAGAAAATCAACAAGATAGACAAGCCCTTATCCAAACTAACCAAAAGACAGAGAGAGAGCATCCAAATCAACAAAATCAGAAATGAAAAGGGGGACATAACAACAGACATTGAGGAAATCCAGAGAATCATCAGGTCATACTTCAAAAACCTCTATTCCACAAAACTGGAAAACCTAAAAGAAATGGATAATTTTCTGGATAGTTACCACATACCTAAATTAAATCAAGACCAGATAAACTATTTAAATAGTCCAATAACCCCTAACGAAATAGAAACAGTCATTAAAAGTCTCCCAACCAAAAAAAGCCCAGGACCAGATGGTTTCAGTGCAGAATTCTACCAGATCTTCAAAGAAGAGTTAATACCAATACTCTCTAAATTGTTCCATACAATAGAAACAGAAGGAACATTACCAAACTCCTTCTATGAGGCTACAATTACCCTGATTCCCAAACCAAACAAGGATACAACAAAGAAAGAGAACTACAGACCGATCTCCCTCATGAACATTGATGCAAAAATACTCAATAAAATACTGGCAAACAGACTCCAAGAACACATCAAAACAATTATCCACCATGATCAAGTAGGATTCATTCCAGGGATGCAAGGATGGTTCAACATACGAAAGTCTGTCAATGTGATACACCATATAAACAAACTCAAAGAAAAAAACCACATGATCATCTCACTAGATGCTGAAAAGGCATTTGACAAAATCCAACACCCCTTCATGATAAAGGTCTTGGAGCGATCAGGAATACAGGGAACATACCTAAACATAATAAAGGCAATTTACAGCAAGCCAACAGCCAACATCAAATTAAATGGAGAGAAACTCAAAGCAATTTCACTAAAATCAGGAACGAGGCAAGGCTGTCCGCTCTCCCCATACTTATTCAATATAGTACTTGAAGTTCTAGCCAGAGCAATAAGACAACATAAGGAGATTAAGGGGATACAAATCGGAAAGGAAGAAGTCAAGCTTTCCCTATTTGCAGATGACATGATAGTATACTTGAGCAACCCCAAAGATTCCACCAAGGAACTGATACAACTTATAAACACCTTCAGCAACATAGCAGGATACAAGATCAACTCAAAAAAATCAGTAGCCCTCCTATATACAATGGACAAAAAAGCGGAGAAGGAAATCAGAGATACATCACCCTTTACTATAGCCACAAATGACATAAAATACCTTGGGGTAACACTAACCAAGCAAGTGAAGGACTTATATGACAAGAACTTTAAGTCCCTGAAAAAAGAAATTGAAGAAGATGTCAGAAAATGGAAAGATCTCCCATGCTCATGGATAGGCAGAACTAACATAGTAAAAATGGCAATCTTACCAAAAGCAATCTACAGATTCAATGCAATCCCCATCAAAATACCAACACAATTCTTCTCAGACCTGGAAAGAATAATACTCAACTTCATATGGAAAAACAAAAAACCCAGGATAGCCAAAAGAATCCTGTACAATAAAACAACCTCTGGAGGCATCACGATCCCTGACTTCAAGCTGTACTATAGAGCTACAGTAATAAAAACAGCTTGGTACTGGCATAAAAACCGACATGTGGACCAATGGAATCGAATTGAAGACCCTGACATTAATCCGCACACCTATGAACAAATAATTTTTGACAAAGAAGCCAAAAGTGCACAATGGAAAAAAGAAAGCATCTTCAACAAATGGTGCTGGCAAAACTGGATATCAACATGTAGAAGGCTGCAAATAGATCCATATCTATCACCGTGCACAAAACTTAAGTCCAAGTGGATCAAGGACCTCAACATAAATCCAGCTACTCTGAACCTGCTAGAAGAGAAAGTTGGAAGTAGTCTTGAACGCATTGGCATAGGAGACCACTTTCTAAATAGAACACCAGTAGCACAGACACTGAGAGAAACAATCAATCAATGGGACCTCTTGAAACTGAGAAGCTTTTGTAGGGCAAAGGATACGGTCAACAAAGCAAAGCGACAGCCTACAGAATGGGAAAAGATATTCACCAATCCCACATCTGACAGAGGACTGATATCCAGAATATATAAGGAACTCAAGAAATTAGACACCAAAATGCCCAACAGTCCAATTAAGAAATGGGCTATAGAACTAAACAGAGAATTCTCAACAGAGGAAACTCAAATGGCTGAAAGACATTTAAGGAATTGCTCAACATCCCTAATCATCAGGGAAATGCAAATCAAAACAACTCTGAGATACCACCTTACGCCTGTCAGAATGGCTAAGATCAAAAACACTGAAGACACCTTATGCTGGAGAGGATGTGGAGCTAGGGGAACTCTCCTCCACTGCTGGTGGGAATGCAAGCTTGTACAACCACTTTGGAAATCAATATGGCGATTTCTTAGAAAATTGGGAATCCATCTCCCCCAAGATCCAGCTATACCACTCTTGGGCATATACCCAAGGAATGCTCAACCACACCACAAGAGCACTTGTTCAGCTATGTTCATATCAGCATTGTTTGTAATAGCCAGAACATGGAAACAACCTAGATGCCCTTCAACTGAAGAATGGATAAACAAAATGTGGTACATATACACAATGGAATACTACTCAGCAGAGAAAAACAATGACATCATGAGGTTTGCAGACAAATGGATGGATCTAGAAAAAATCATCCTGAGTGAGGTATCCCAGACTCAGAAAGACAAACATGGTATGTACTCACTCATAACAGGATACTAGATGTGGAACAAGGATGACTGGACTGCTACTCACATCACCAGGCAGGCTACCTGGAAAACAGGACCCCAAGAAAGACACAGGGATCGCCCAATGACAGAGAAATGGAATGAGATCTACATGAACAGCCTGGACATGAGTGGGGGTAGTGAAGGGCGAAGGTCGAGGGAAAGAGAGCTTGGGTGAGTGGGAGATCCCAGCTGGATCAACAACAGAGAGGGAGAACAAGGAATAGGAGACCATGGTAAATGAAGACCACATGAGAATAGGAAGAAACAAAGTGCTAGAGAGGCCCACAGAAATCCACAAAGATACCCCCACAACAGACTGCTGGCAATGGTCGAGAGACAGTCCGAACTGACCTACTCTGGTGATGGGATGGCCAAACACCCTAATTGTCGTGCTAGAAACCTCATCCAACTACTGAGGGATCTGGATGCAGAGATCCATGACTAGGCCCCAGGTGGATCTCTGGGAGTCCAATTAGCGAGAATGAGGAGGGTTTATATGAGCGAGAATTGTTGAGACCAAGGTCGGATTAAGCACAGAGACAAATAGCCGAACGAACGGAAACACATGAAATATGAACCAATGGCTGAGGGGTCACCAACTGGATCAGGCCCTCTGAGTGGGTGAGACAGTTGATTGGCCTGATCTGTTTGGGAGGCATCCAGGCAGTGGCACCGGGTCCTGTGCTCATTGCATGAGTCGGCTGTTTGAAACCTGGGGCCTATGCAGGGTCCCTTGGCTCGGCCTGGGAGGAGGGGACTGGACCTACCTGGACTGAGTCCACCAGGTTGATCTCAGTCTGTGGGGAAGGCTTTGCCCTGGAGGAGATTGGAATCGGGGGCGGGCTGGGGGGGAAGGTGAGGGAGGCGGGAGGGGGGAGAACAAGGGAATCTGTGGCTGATATGTAGAACAAGGGAATCTGTGGCAGATC

I did a BLASTn of this sequence against the Mus database to confirm that this sequence (from my testes
transcriptome) does in fact match to the gene Lhcgr

The BLASTn result was NOT Lhgcr:

Peromyscus maniculatus beta-globin gene cluster, complete sequence	8576	21368	84%	0.0	91%	EU204642.1
Select seq gb|EU559333.1|	Peromyscus leucopus beta-globin gene cluster, partial sequence


***THEREFORE I WILL NOT BE USING THIS TRANSCRIPT ID TO SEARCH FOR LH IN THIS ANALYSIS***


---

***Lhb serach***

grep "BINPACKER.9392.2" good.BINPACKER.cdhit.fasta -A 1


>BINPACKER.9392.2
CAGAGAGTGTGAGAACTATGGCAGGCCGGACACTAGCTCTGCGATATGGACCCCCTTGGTCCCCCATTTCTGAAACTGAGGTTCCTGGAACTTGGCCCAGCTGGCATCTTACCAGCAGCGGGGTCGCCCACCACAGGATCCCACCGGCACCCTTTCCTCCGCCCATTCTACAGCCTACGGTCCTGGCACCCCTGCCCGCAGCCGTGAGGCACGACCCGCGCATCTGGGCCTTTGACGAAGTCATCAACAGATGGGAAACCACCTCAGGCTCGGCACACAGGCCCAAGACACACAGCGGGCCCTGTGCTCAGCCCAAGGCAGCAGAGCATGAGGACCCTGGGCGGGTACTTGGGATCAAGTCTTTGGCAGACAAGCTAAGAAGACACCAGGGTTGGGGTGTCCCCCTGGACACAAAATACCAGATCAGCGAGACGAAGGCCGAGTACACGGGCTGCCCAGGTCTGGAGCAGAGTGCCCCTCTCTTTGTGGAGCCCCAACCCCCGGAGCTCGCTGACCACCACCGTGGGGGCCCATCCCAGGCCCTAATCCCCTGGACGAGAAACCCAGAGCTGGCTGGCCAGCCCTTTACAGTGTGTAAGATGGGTGTCCTGGGCCGCCTCCAGCCCTATCTGACCACCTCCGTACGTGACTTCTCCAGGTGAGCCTCTGGCCCTCTCCTCTGCAGTTCCAACTGCAGATATAGGAACCGCAGTACTGACCCAAGCTCCACTCTGGAGCTCAGAACCCAGAGCCTCCTCTAAGGCAGCAAGAGCCCATTTCCAGGATCCGAACTCTGTTGCTCTTGGTCCTATTCACCTGCTGTTCCACCCACAGGAAAGAGCTGTCTGGACACCCTGGCCGGGACACAGTGGTGAAGTCACAGCGCCCCCGCCACCTTAAGCGGCTACCAAGGGAGCGCCTAGCTCGGGCGCGCCCAGTGCCACCGGCTGTGCCATACCTAGGGGCCCTGCCTCTGACTCAGGAGTCTTATGGGCCCTCGGTGCACCCGCTCCGCAAGCTGGACCGCTTTTGCCCACTAGAGGCCCCCTGGGGAGGCCCCCACCCGAAGCCAGTGCCTGGCATCTACAGCGTGCCCAAGGCCTACTGCACTGAGAACTCCCGCTATGGGAGTGCCAGGGCAGAGCTGGTGTGAGCCCACCAGGCAGACATACACTCCAGCATCGGATTTGTGGGTTGTGGGCAGGTCAACCCAGGCATACCCCATGCTTGCAGAGGAGCCTCCTCCAGGCGGAGCTTACCAATCAAGGAGTTTTATAGCTGAAACCACACCCATTTCTGGACCCATGCATCCTGATTAGGGAGTTGGGTGAGGGGTGGTTGCTCCACCTCTGGTTGAAGGTGTCTGGGTATAGTTGGAGGCCCACTTGGCATGGGGGGAACCTAAGGCCATAGCTGTACCTTTTATAAGGTTGCACCGAGGCATGAGAATGAGCTAGCCTCCAGCACCTTAGGGGCCCAGCCAGCAGCCTACAGGATACTCTCCCCTTTCCCCTTCCGTGATGGGGTAAAGGATCCAAGGGTCCTTCCTGCACTCCCAATGTCCGCTAGGCTCACACCTGGGCTGAGTATGAGGCCGATCACCGTGACACAGGAGCTGGTCCCTCACTTCCCTGACCTTGTCTGTCACCGCCCCCAAAGAGATTAGTGTCTAGGTTACCCAAGCCTACAGCCGCTGCTTGGTGGCCTTGCCACCCCCACAACCTGCAGGTATAAAGTCAGGTGGCCAAGGTAGGGAAGGCATCGAGAATGGAGAGACTCCAGGTAAGAGTGTAGGGTCCGGGACATCTCCCAACTCCACCAGACCCTGGCTGGAACAGATGGACAGCCCTGTGACCTGGGGGGAGGGGCGTCCTAGCTGGTGGGCTCCCAGCAGGACAGAAATGGATTGGATGGCAGGTGATGGGTCCTGAAGGCCAGGATCTATATTACACTGAATGGGTCCAGAGCCCGGGATGGAAACCCAGGGCTAGACTGAGACACTGGCTATGTCCCAGGGGCTGCTACTGTGGCTGTTGCTGAGCCCAAGTGTGGTATGGGCCTCCAGGGGCCCCCTGCGGCCACTGTGCCGGCCTGTCAACGCAACCTTGGCCGCAGAGAATGAGGTCTGCCCAGTCTGCGTCACCTTCACCACCAGCATCTGTGCCGGCTACTGTCCTAGCATGGTTCGAGTACTGCCGACTGCCTTGCCTCCTGTGCCCCAGCCTGTGTGCACCTACCGTGAGCTGCGCTTCGCGTCTGTCCGCCTCCCTGGCTGCCCACCTGGTGTGGACCCCATGGTCTCCTTCCCCGTGGCCCTCAGCTGCCGCTGTGGGCCCTGCCGTCTCAGCAGCTCTGACTGTGGGGGTCCCAGGGCTCAACCGATGGCCTGTGACCTCCCCCACCTCCCCGGCCTCCTCTTCCTCTGAGGCCCACCCGCTAACTCCTCTATCTCAAGGAGCTGGCACGTGTTCCAACTGTCCCTCCCAATAAAGGCTAATTTACAACTGC

I did a BLASTn of this sequence against the Mus database to confirm that this sequence (from my testes
transcriptome) does in fact match to the gene Lhb

The search results confirm this is a good Lhb sequence match:

Mus musculus adult male hypothalamus cDNA, RIKEN full-length enriched library, clone:A230101A12 product:luteinizing hormone beta, full insert sequence
Select seq ref|XM_006228975.2|	PREDICTED: Rattus norvegicus luteinizing hormone beta polypeptide (Lhb), transcript variant X2, mRNA	1411	2103	89%	0.0	79%	XM_006228975.2
Select seq gb|U25653.1|RNU25653	Rattus norvegicus testicular luteinizing hormone beta-subunit (TLHB1) mRNA, complete cds


***Therefore the Lhb sequence in Mus corresponds to the seq >BINPACKER.9392.2 in the Testes transcriptome
AND I WILL USE THIS TESTES TRANSCRIPT ID FOR THIS SEQUENCE FOR MY LH ANALYSIS***


***I will see if this transcript (>BINPACKER.9392.2) corresponds to a gene in my gene x transcript matrix***  

Note: NEWESTFinalMUS.txt is file that only has one transcriptID per MUS gene ID that has matched

The matrix file is located here:

lauren@davinci:/mnt/data3/lauren/TESTESdge/MUS/NEWESTFinalMUS.txt

grep "BINPACKER.9392.2" NEWESTFinalMUS.txt

RESULT:

BINPACKER.9392.2	ENSMUSG00000100916.3


This match is reflective of the following ensemble entry:


ENSMUSG00000100916 luteinizing hormone beta (Lhb)


So I can search for this "ENSMUSG00000100916.3" entry in my dataset within the EdgeR DGE analysis



***I graphed this result in with the other nine genes I previously graphed***


tag.keep<-c("ENSMUSG00000100916.3", "ENSMUSG00000079019.2","ENSMUSG00000001768.15", 			"ENSMUSG00000054200.6",
            "ENSMUSG00000026435.15","ENSMUSG00000020427.11","ENSMUSG00000019997.11",
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




	d. MERGE QUANT FILES
	
***The purpose of this merging step is so that the data is compiled properly for DTE***
	
All of the quant files are in one folder, and I used the following MacManes script to merge the quant files into one:

The merge_quants.py script was saved in the salmon folder, where all of the quant outputs are.


* This file is in the following folder:

lauren@davinci:/mnt/data3/lauren/TESTESdge/NEWmaps/salmon/merge_quants.py


The following is the merge_quants.py script:

!/usr/bin/env python

USAGE: python merge_quants.py dry/ dry.counts

--
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

--

The command to run the script was:

python merge_quants.py /mnt/data3/lauren/TESTESdge/NEWmaps/salmon/ NEWmerged.counts

This file was also downloaded and named NEWmerged_counts.csv


The data appear with the transcript ID as the left most column, and then each sample is an additional column, with the raw transcript counts within each sample for that transcript ID.  
Thus each row contains the unique transcript ID with the counts of transcripts (raw counts, NOT TPMs) for each sample.  



	e. differential transcript expression with edgeR


***DEPENDENCIES***

source("https://bioconductor.org/biocLite.R")

biocLite("tximport")

library(tximport)

install.packages("readr")

library(readr)

biocLite("edgeR")

library(edgeR)


***MY DATASET:***

***My dataset is NEWmerged.counts (and it is csv form) and I import it below:***


x <- read.csv("~/Desktop/NEWSALMONdge/NEWmerged.counts", header=TRUE, row.names=1)
head(x)


* only 1 sample needs a cpm >1 to keep the transcript ID , this  seems most logical:

z <- x
keep <- rowSums(cpm(z)>1) >= 1
y <- z[keep,]
dim(y)
y[y==0]<-1


***my groups:***

group <- factor(c(1,1,1,1,2,2,2,2,1,2,1,1,1,1,2,1,2,2,2,1,2,2))


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


***These are my results:***

   [,1]
    
-1    45

0  69021

1     21

* 66 total significant genes
* 45 higher in WET
* 21 higher in DRY


***These are the BinPacker transcript IDs for the significant genes:***

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



***This is a command to print out the 66 significant results:***

topTags(et, n=66, adjust.method="BH", sort.by="PValue")

***RESULTS:***

Comparison of groups:  2-1 

                       logFC      logCPM       PValue          FDR
                       
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


***plots***


plotSmear(et, de.tags=detags, main="DE genes, all data, FDR 0.05")
abline(h = c(-2, 2), col = "blue")
abline(h = c(-4, 4), col = "blue")

-----------------------------------------

***Now that I have the transcript IDs which are significantly differentially expressed between groups, 
the following will allow me to quickly view side by side the transcript ID, gene ID, & gene info:***

cat NEWESTmusHITS.txt | awk '{print $1 "\t" $7 "\t" $10 "\t" $11}' > NEWESTGENEmusHITS.txt

cat NEWESTGENEmusHITS.txt | wc -l

* 45636

cat NEWESTGENEmusHITS.txt | sort -uk1,1 > NEWESTSORTEDgeneHITS.txt

cat NEWESTSORTEDgeneHITS.txt | wc -l

* 37744


***Then I can view the gene ID and gene info for each transcript ID that was significant:***

grep "BINPACKER.9726.2" SORTEDgeneHITS.txt


***I also retrieved the BINPACKER transcripts which did not have a matching Ensembl ID in the DTE:***

This was done by the following command in the folder that contains my transcriptome assembly:

/mnt/data3/lauren/TESTESdge

example: 

grep "BINPACKER.3452.1" good.BINPACKER.cdhit.fasta -A 1

>BINPACKER.3452.1
CTCCAAGGGATCGGCTTCAAATGGACCTATAACACAGTGATGCGTCCAGGAGAACCAGTGCCCATCTCTAACACTAACAAAGTGAAGGAAGATATCAAATACAAATACAAAGTTGAAGACGATGATGATGATGAGGATGATGATGATGATGAGGATGATCATAACATTGATGAGGAGTAGAAAGATGAGAAAAGGTAGTAGAGAGCTCAGATCTCAGTCCTGTATCGCCACATTGTATGACCTCCTCTCAGCCAAGCTCACCCATCTGGAGGGACCAGTGTACCAGGTGACCCTGGGCTGGAGGAGGCAAAGACCAAATGGCCCTATATCTCTGTGTCCAATAAACTCTGTGCATCCATATGAGCCTGAGTGGTCCTTGTCACCAAATGGGGTGGTTTGAGGTCCCCCAGAGCTCCTGTACCTCTTCACTTTTTGGACCTACACCCGTGCACTGAGTTAGGGTTTACACAGAAAGACCCTTGGTCAACTCAGCAGGCCTGAGGTCCCCAGAACGAGGCCAGCATCCAGCCCTGCCTGTCTTAAAGCCACACAGGGCCTGCAGGAAGAGCTCAATGTGTGTTGATCAAGCTGGCGGGGCCTGGCTTTGCAGTGGTCCCATGTAACTGTATTGAGCAGGTGCCATTGTCTCCTGTTCCTCTCCAGACCCTCCCGTTCTCCAGCCACCTCTTCCCTTGCCTTCTCCTGTCCTGGTCACTATCCTCCCTCACAAAGCTAGGCTGCCTGAAATGCATCTCCTAAGCAGTTGGGGATTAAATCTAAGGCCTCATGAATGCTCGGCAAGTACTCTAGCACTGAACCGTATTCTTAGCTTGCATTTATTTTTTTAGTATTCATGAATAAATGAATGAATGAATGATACTTATATATATGCAAGCAAACCACTCGCACACCTA

* All if these grepped sequences and their results are in a different file: 
***DTEno-matchBLASTnSequences.md***

	e. differential gene expression with tximport and DESeq2


***DEPENDENCIES:***

source("https://bioconductor.org/biocLite.R")

biocLite("tximport")

library(tximport)

install.packages("readr")

library(readr)

biocLite("DESeq2")

library(DESeq2)



***import the gene ID matrix***

tx2gene <- read.delim("~/Desktop/NEWSALMONdge/NEWESTFinalMUS.csv", header=FALSE)
head(tx2gene)


***SET DIRECTORY***

dir <- setwd("/Users/lauren/Desktop/NEWSALMONdge/test")
getwd()
file.path(dir, "salmon", "quant.sf")
dir <- setwd("/Users/lauren/Desktop/NEWSALMONdge/test")
getwd()
file.path(dir, "salmon", "quant.sf")

* This is what should be returned by these commands:
      "/Users/lauren/Desktop/NEWSALMONdge/test/salmon/quant.sf"  
      

***SET samples in downloaded order***

samp <- c("salmon_13T", "salmon_102T", "salmon_209T", "salmon_265T", "salmon_343T", "salmon_344T", "salmon_349T", "salmon_355T", "salmon_366T", "salmon_376T", "salmon_381T", "salmon_382T", "salmon_383T", "salmon_384T", "salmon_400T", "salmon_888T","salmon_999T", "salmon_1357T", "salmon_1358T", "salmon_1359T", "salmon_2322T", "salmon_3333T") 
samp


file.path(dir, "salmon", samp, "quant.sf")
files <- file.path(dir, "salmon", samp, "quant.sf")
files

names(files) <- samp
names(files)


***NOW read in files with tximport***

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, read = read_tsv)
head(txi.salmon$counts)


***Read in table for group names

sampleTable <- read.csv("~/Desktop/mergedCounts/sampleTable.csv", header=TRUE, row.names=1)
head(sampleTable)

* Note: Dry = 2 ; Wet = 1

sampleTable$condition <- as.factor(sampleTable$condition) 

* (assigned the condition column as a factor)

rownames(sampleTable)==colnames(txi.salmon$counts)

* (check if all true, and all are true, so proceed)



dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~condition)

* (prepare the entire matrix for analysis by condition)


dds <- dds[ rowSums(counts(dds)) > 1, ]

* (filter out the rows-transcript IDs-with less than 1 read) 

* Also note that the count data is counts(dds)

dds <- DESeq(dds)

* (perform DEseq statistical analysis)



res <- results(dds, alpha=0.05, pAdjustMethod = "BH")

* (results filtering based on normalized counts of each gene with a FDR value set at alpha= 0.05 for adjusted p value)

summary(res)


***RESULTS:***

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

* 215 significant results


resOrdered <- res[order(res$padj),]

* Order the results from most to lease significant by adjusted p-value


write.csv(resOrdered[1:total,],file="NEWresults.csv")

* All of the results were written to a file so that I could view them


plotMA(resMLE, MLE=TRUE, main="unshrunken LFC", ylim=c(-4,4), colLine=NONE)
abline(h = c(-2, 2), col = "blue")
abline(h = c(-4, 4), col = "blue")



	f. Comparison of Log2FC values for DGE in EdgeR with DGE in SESeq2
	

***Extract log2 FC values for each ensembleID in the et matrix from EdgeR***


output <- topTags(et, n=14218)

write.table(output, file="diffexp_detags_edgeR.csv", sep = "," , row.names = TRUE)


* within et (my table), n = the total number of ensembl gene IDs


***Extract log2 FC values for each ensembleID in the res matrix from DESeq2***


write.table(res, file="diffexp_DESEQ2.csv", sep = "," , row.names = TRUE)

* my data table is res



***NEXT Extract data frame from edgeR result object (from the list which is "output" select the .Data item in the list and take the first element of .Data, which is the results table, and convert this to a dataframe, named edgeRtable:***

edgeRtable<-data.frame(output@.Data[1])


***Add rownames (Gene ID) as a column in each table***

edgeRtable$ID<-rownames(edgeRtable)

res$ID<-rownames(res)


***Keep only the columns wanted and convert Fdeseq to dataframe***

Fdeseq<-as.data.frame(res[,c("log2FoldChange","ID")])

FedgeR<-edgeRtable[,c("logFC","ID")]



***Merge final edgeR table with final DESEQ table***

compare<-merge(Fdeseq,FedgeR,by="ID")


***Change the column names so that we label each lFC by the analysis***

names(compare)<-c("ID","D_l2FC","E_l2FC")


***linear regression of results (so that we can plot with a regression line)***

regression<-lm(compare$E_l2FC~compare$D_l2FC)

summary(regression)


	Call:
	lm(formula = compare$E_l2FC ~ compare$D_l2FC)

	Residuals:
     Min       1Q   Median       3Q      Max 
	-2.38746 -0.09001 -0.00879  0.09039  2.32470 

	Coefficients:
               Estimate Std. Error t value Pr(>|t|)    
	(Intercept)    0.305650   0.001669   183.1   <2e-16 ***
	compare$D_l2FC 1.380585   0.008319   166.0   <2e-16 ***
	---
	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

	Residual standard error: 0.1974 on 14214 degrees of freedom
	Multiple R-squared:  0.6596,	Adjusted R-squared:  0.6596 
	F-statistic: 2.754e+04 on 1 and 14214 DF,  p-value: < 2.2e-16



***Plot with regression line and axis labels (regression line optional)***

plot(compare$D_l2FC,compare$E_l2FC,ylab="EdgeR_lfc",xlab="DESeq2_lfc")
abline(regression)


***Correlation of results***

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




