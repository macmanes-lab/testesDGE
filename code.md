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

n	n:500	L50	min	N80	N50	N20	E-size	max	sum	n...
