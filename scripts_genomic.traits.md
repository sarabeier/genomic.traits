# Scripts and comand lines for the computation of genomic traits and species-level taxonomy

### by Sara Beier with support of Johannes Werner

This file contains unix shell command lines and R scripts to compute genomic traits as well as the species-level taxonomic affiliation based on the GTDB taxonomy or based on the average nucleotide identity
from the genome data as given in Table S1 in the manuscript Beier et al, doi: https://doi.org/10.3389/fmicb.2022.985216. The information in the folloing columns of Table S1 is provided via the JGI genome statistics (https://img.jgi.doe.gov/), : NCBI_genus, Cultured, Habitat, NCBI.Assembly.Accession,NCBI.Bioproject.Accession, NCBI.Biosample.Accession, NCBI.Taxon.ID, Genome.size, Gene.count, %GC, RRN_IMG, %HGT. Information in the remaining columns was obtained as detailed below.




## rrnDB 16s rRNA copy number (RRN_rrnDB)

Analyses indicated that 16S rRNA gene copy numbers (RRN) available via the JGI/IMG platform may in many cases underestimate true values for RRN (see Figure S1). For our trait covariation analyses we have therefore extrapolated RRN values from the rrnDB ([rrnDB](https://rrndb.umms.med.umich.edu/), Stoddard et al. 2015) , a database that provides manually curated RRN values. Genus level RRN values from individual rrnDB entries (mean values ) were computed using the Rscript below and matched with the genus information (NCBI_genus) given for each entry of the JGI/IMG genomes provided in Table S1.

```R
#rrnDB data
rrnDB <- read.csv ("rrnDB-5.7.tsv", header=T, sep='\t', fill=T)[,c(1:4,7,12)] #dataset available via https://rrndb.umms.med.umich.edu/static/download/
rrnDB.genus <- tibble(rrnDB) %>%
    filter(grepl("(genus)",RDP.taxa)) %>% #selects entries with genus level affiliation
    mutate_if(is.character, str_replace_all, pattern = ' \\(genus\\)', replacement = '%') %>% #replaces ' (genus)' with '%'
    separate(RDP.taxa,c("genus","rest"), sep="%") %>%  #splits column by separator % 
    select(-c(1,2,3,4,6)) %>% 
    group_by(genus) %>% 
    summarise_at(vars(X16S.gene.count), mean, na.rm = TRUE) 
```
Citation:  
* Stoddard, S. F., Smith, B. J., Hein, R., Roller, B. R. K., and Schmidt, T. M. (2015). rrnDB: improved tools for interpreting rRNA gene abundance in bacteria and archaea and a new foundation for future development. Nucleic Acids Res. 43, D593–D598. doi:10.1093/nar/gku1201.


## Read-based estimations of RRN and genome size
RRN estimates given in the rrnDB as well as in the JGI/IMG platform are obtained from assembled genomes. The detection of RRN after genome assembly may however suffer from underestimations, because the 16s rRNA gene is highly conserved and multiple copies that are (nearly) identical may collapse into a single copy during assembly. This biases should be less important if estimating the full genome size (and other traits) from assembled data, as the 16s rRNA gene typically represents <0.5% of the whole genome and protein-coding genes are usually less concerved and accordingly less sensitive to assembly errors.

An alternative assembly independent option to estimate RRN counts is to relate the number of reads encoding the 16s rRNA gene to the number of reads coding for obligatory single copy genes in a genome (Biers et al., 2009). Similarely, the genome size can be estimated by relating the total number of reads to the number of reads encoding obligatory signle copy genes in a genome (Nayfach and Pollard, 2015). 

We could not extract sufficient original read data for genomes used in this study to relay the main analyses of this study on read based RRN and genome size estimates. However, we could compare assembly based RRN and genome size values from the JGI/IMG or the rrnDB databases against read-based estimates for ~1500 genomes. While also the read-based estimation of RRN or genome size may suffer from biases, these biases are not assembly-dependent and should accordingly not covary with biases of RRN estimates after assembly. We therefore interpreted the strength of correlation between assembly-based and read-based estimates for RRN and genome size as measure for the quality of values provided in the rrnDB and JGI/IMG databases.

### Extraction of read accession numbers
This script extracts SRR runids from the NCBI Bioproject Accession available in Table S1. Duplicate NCBI Bioproject Accession numbers in Table S1 were excluded as here several genomes may have been sequenced here under a single NCBI Bioproject Accession. Reads from metagenome assembled genomes were excluded because reads from a single genome cannot be seperated out.

```bash
sort TableS1.tsv |grep -v 'Metagenome-Assembled Genome'|cut -f15 |uniq -u| while read line 
do
wget -qO- "https://www.ebi.ac.uk/ena/browser/api/xml/${line}" | grep -A1 "ENA-RUN" | sed '1d' | sed -e 's/^ *<ID>//' -e 's/<\/ID>//' | paste <(echo ${line}) -
done >> bioproject_runid.txt
```
The extracted SRR runids were used to download raw read files. To speed up the downstream analyses we used only the reverse run file for paired-end reads and all read files were subsamples to maximal 2.000.000 reads. Read files ending with *_consensus.fastq.gz or _subreads.fastq.gz were excluded, as these files result from sequencing methods producing long reads that may include more than one gene.

### Determination of reads encoding the 16s rRNA gene
We used the SortMeRNA software (Kopylova et al., 2012) to identify reads encoding the 16s rRNA gene. Subsequently we counted the number of 16s encoding reads, all reads and determined the average read length within one sample.

```bash
sortmerna \
--workdir . \
-ref silva-bac-16s-id90.fasta \
-ref silva-arc-16s-id95.fasta \
-reads sub.SRRrunid.fastq.gz \
--fastx true \
--other false \
--threads 1 \
--no-best true \
--num_alignments 1

gunzip out/aligned.fq.gz
echo $(cat out/aligned.fq | wc -l)/4|bc > 16s.reads #number of total reads
echo $(zcat ../SRRrunid.fastq.gz | wc -l)/4|bc > tot.reads #number of reads encoding the 16s rRNA gene
echo SRRrunid > names #assigns runids to the output data
cat out/aligned.fq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > aligned.fasta
perl -e '$count=0; $len=0; while(<>) {s/\r?\n//; s/\t/ /g; if (s/^>//) { if ($. != 1) {print "\n"} s/ |$/\t/; $count++; $_ .= "\t";} else {s/ //g; $len += length($_)} print $_;} print "\n"; warn "\nConverted $count FASTA records in $. lines to tabular format\nTotal sequence length: $len\n\n";' aligned.fasta  > aligned.tab
perl -e ' $col=2; while (<>) { s/\r?\n//; @F = split /\t/, $_; $len = length($F[$col]); print "$_\t$len\n" } warn "\nAdded column with length of column $col for $. lines.\n\n"; ' aligned.tab > l.aligned.tab
awk -F'\t' '{sum+=$4} END {print "l.aligned.tab", sum/NR}' l.aligned.tab > read.lengths.16s #returns average read lengths
```
if looping the commands above through all read files, the results can be summarized with the following command:
```bash
paste names tot.reads 16s.reads read.lengths.16s >> res.tab
```
### Determination of read-based genome size and the number of sequenced genome equivalents
We used the software MicrobeCensus (Nayfach and Pollard, 2015) to estimate genome size from read data of sequenced genomes applying the following command:
```bash
run_microbe_census.py -t 1 SRRrunid.fastq.gz genome.size.SRRrunid
```
After looping through all read input files the results can be summaried using the following code

```bash

grep 'average_genome_size:' genome.size* >genomesize.sum.tab
grep 'genome_equivalents:' genome.size* >genomeequivalents.sum.tab

paste genomesize.sum.tab genomeequivalents.sum.tab |sed 's/_size:    \'$'\t''/_size:\'$'\t''NA\'$'\t''/g'|sed 's/_equivalents:    /_equivalents:\'$'\t''NA/g' |sed 's/genome.size.sub.//g'|sed 's/:average_genome_size://g' |cut -f1,2,4 | sed 's/\.R/\'$'\t''R/g'>mc.out.tab
```
A common summary result file from the SortMeRNA and MicrobeCensus output was created to estimate read-based RRN:

```bash
paste mc.out.tab res.16s.tab |cut -f 1-4,7-9 > genomesize.16s.tab
```
### Determination of read-based RRN
Read-based RRN was estimated in R by normalizing the number of reads encoding the 16s rRNA gene by read lengths relative to the lengths of the total 16s rRNA gene (in E. coli: 1541 base pairs) and divide the resulting value through the number of sequenced genome equivalents:

```R
# Load libraries
library(tidyverse)

# Load input data
l16s <-1541 #full length of E. coli 16s rRNA gene
gs16s <- read.table("genomesize.16s.tab", sep = '\t') #combined SortMeRNA and MicrobeCensus outputfile
rownames(gs16s) <- paste(gs16s$V1, gs16s$V2, sep='.')
colnames(gs16s) <- c('runID','readID','genome.size.MC','genome.eq.MC','n.reads','n.16s','l.16s')
gs16s[gs16s == 0] <- NA
IDs <- read.table("/Users/sara/Documents/R-scripts/GenomicTraits/data/bioproject_runid.txt",  sep = '\t')
colnames(IDs) <- c("NCBI.Bioproject.Accession" ,"runID")
gtraits <- read.table('/Users/sara/Documents/DFG/GenomicTraits/JGItraits2021/TableS1.tsv',  sep = '\t', header=T)


dat <- tibble(gs16s) %>%
    inner_join(IDs, by = c("runID" = "runID")) %>% #merge with IDs, keep only row in both dataframes
    drop_na %>% #remone rows with NA
    filter(n.reads>=2000000) %>% #filters out all genomes with < 2.000.000 reads
    filter(genome.eq.MC>50) %>% #filters out all genomes with < 50 genome equivalents sequenced
    mutate(RRN_read.based = round(n.16s *l.16s/l16s/genome.eq.MC,0)) %>% #number or reads normalized per read length and tot 16s length divided by genome equvivalents
    filter(RRN_read.based<30) %>% #for 9 genomes a readbased RRN > 30 was detected. We considered this as not reasonable values and removed this 9 genomes
    inner_join(gtraits[,c(1,9:12,20,21,24,25,26,27,29,31)], by = c("NCBI.Bioproject.Accession" ="NCBI.Bioproject.Accession")) %>%
    mutate(taxon_oid=as.character(taxon_oid)) %>%
    filter(runID != 'SRR1810285') #strain with one single runID but multiple assemblies >> multiple picrustIDs with different RRNs, genome removed

# Create Table S2
TableS2 <- dat[,c(13,1,22,24,11,19,3:4,17)]
colnames(TableS2) <- c('IMG.Genome.ID','runID','RRN_IMG','RRN_rrnDB','RRN_read.based',
'Genome.size', 'Genome.size_read.based', 'Genome.equivalents_read.based','GTDB_genus')
```

Citations:  
* Kopylova, E., Noe, L., and Touzet, H. (2012). SortMeRNA: fast and accurate filtering of ribosomal RNAs in metatranscriptomic data. Bioinformatics 28, 3211–3217. doi: 10.1093/bioinformatics/bts611.
* Nayfach, S., and Pollard, K. S. (2015). Average genome size estimation improves comparative metagenomics and sheds light on the functional ecology of the human microbiome. Genome Biol. 16, 51. doi: 10.1186/s13059-015-0611-7.
* Biers, E. J., Sun, S. L., and Howard, E. C. (2009). Prokaryotic genomes and diversity in surface ocean waters:interrogating the Global Ocean Sampling metagenome. Appl. Environ. Microbiol. 75, 2221–2229.


## Genome data accession

The remaining genomic traits were computed based of sequence or annotation data provided for genomes via the the JGI/IMG platform, which can be downloaded as tar-compressed files. The scripts below demontsrate trait computations for the example genome 2593339298 (Mameliella alba DSM 26384, [JGI/IMG platform](https://img.jgi.doe.gov/), 2593339298.tar.gz)

```bash
$ tar -xzvf 2593339298.tar.gz

$ ls 2593339298
2593339298.cog.tab.txt  2593339298.genes.fna       2593339298.ipr.tab.txt   2593339298.signalp.tab.txt  README.txt
2593339298.fna          2593339298.gff             2593339298.ko.tab.txt    2593339298.tigrfam.tab.txt
2593339298.genes.faa    2593339298.intergenic.fna  2593339298.pfam.tab.txt  2593339298.tmhmm.tab.txt
```
## Genome richness and gene duplication
This code returns the number of different gene orthologs within a genome (gene richness) and the averagy copy number of all gene orthologs within a genome (gene duplication). 
Input file: JGI COG annotation table
```R
cog <- read.table(file="2593339298.cog.tab.txt", header =T, sep='\t', fill=T, stringsAsFactors = FALSE, quote = "")[,1:10]
rich.cog <- length(table(cog$cog_id)[table(cog$cog_id)>0]) # Gene richness
dup.cog <- mean(table(cog$cog_id)) # Gene duplication

#rich.cog
#1799
#dup.cog
#1.926626
```

## Transcription factors (%TF)

This code returns total number of transcription factors. To obtain fraction of transcription factors (column %TF in Table S1) divide by the number of total genes (column 'Gene count' in Table S1) given via the JGI platform or in Table S1.
Input file: DBD-PFAM.list, list pfams with DNA binding sites from the DBD database (Wilson et al., 2008).

```bash
$ head DBD-PFAM.list 
#http://www.transcriptionfactor.org/Download/pfam_18.dbds.v2.03.txt
# generated by ./pfamdbddetails.pl 18 pfam_catfile_18
pfam00010
pfam00046
pfam00096
pfam00105
pfam00126
pfam00157
pfam00165
pfam00170

$ TF_count=$(grep -f DBD-PFAM.list 2593339298/2593339298.pfam.tab.txt |wc -l) #count transcription factors
$ echo $TF_count
182
```

Citation:  
* Wilson, D., Charoensawan, V., Kummerfeld, S.K., and Teichmann, S.A. (2008) DBDtaxonomically broad transcription factor predictions: new content and functionality. Nucleic Acids Res 36: D88–D92.



## Prophages

This code returns the total number of category 1 (sure) and category 2 (somewhat sure) prophages detected via the VirSorter software (Roux et al., 2015) (see column'Prophages' in Table S1).

```bash
#run virsorter
$ wrapper_phage_contigs_sorter_iPlant.pl -f 2593339298/2593339298.fna --db 1 --wdir virout.2593339298 --ncpu 4 --data-dir /
virsorter-data/

#extract prophage counts from VirSorter outputfile
$ var1=$(sed -n '/## 4/,$p' virout.2593339298/VIRSorter_global-phage-signal.csv | wc -l) #counts number of lines including and below line: ## 4 - Prophages - category 1 (sure)
$ var2=$(sed -n '/## 5/,$p' virout.2593339298/VIRSorter_global-phage-signal.csv | wc -l) #counts number of lines including and below line: ## 5 - Prophages - category 2 (somewhat sure)
$ var3=$(sed -n '/## 6/,$p' virout.2593339298/VIRSorter_global-phage-signal.csv | wc -l) #counts number of lines including and below line: ## 6 - Prophages - category 3 (not so sure)
$ prophages=$(echo "($var1-$var2-2) + ($var2-$var3-2)"|bc) #returns number of detected prophages from catagory 1 and category 2
$ echo $prophages
4
```

Citation:  
* Roux, S., Enault, F., Hurwitz, B.L., and Sullivan, M.B. (2015) VirSorter: mining viral signal from microbial genomic data. Peerj 3: e985.

## Codon usage bias (CUB(F)) and Generation time (Vieira-Silva)

This code returns the codon usage bias parameter F that was used in the Figures 3 and 4 for CUB, as well as the predicted generation time d.vs as detailed elsewhere (Vieira-Silva 2010)(see column 'Generation time (Vieira-Silva)' in Table S1). The parameter F is based on the earlier defined parameters &Delta;ENC' (Rocha 2004) and S (Sharp et al. 2005). Codon usage bias summary statistics were calculated in a shell script using the ENCprime software (https://github.com/jnovembre/ENCprime). Sequences coding for ribosomal proteins were extracted using the script extractFromFasta.pl (https://github.com/jonbra/NGS-Abel/blob/master/scripts/extractFromFasta.pl). Perl code available via the FAS center page (http://archive.sysbio.harvard.edu/CSB/resources/computational/scriptome/UNIX/) was used to reformat sequence files.
The ENCprime output together with extracted ribosomal gene coding sequences were used as an input for an downstream R script, to return the variables F and d.vs.

```bash
#select annotated genes (without stopp-codon)
$ cd 2593339298
$ awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < 2593339298.genes.fna >enc.outfile.fas
$ awk '{if ( $1 != "") print $1}' <enc.outfile.fas> enc.tmp
$ grep '>' 2593339298.genes.faa |sed 's/>//'g |sed 's/\s.*$//' > enc.2593339298.genes.list
$ perl ../extractFromFasta.pl enc.tmp list enc.2593339298.genes.list  > enc.tmp2
$ sed 's/TAA$//' enc.tmp2 |sed 's/TGA$//'|sed 's/TAG$//'|sed 's/^ATG//' >enc.2593339298.genesnostop.fna #remove final stop-codons and initial start codons

#remove sequences which are not multiple of 3
$ perl -e '$count=0; $len=0; while(<>) {s/\r?\n//; s/\t/ /g; if (s/^>//) { if ($. != 1) {print "\n"} s/ |$/\t/; $count++; $_ .= "\t";} else {s/ //g; $len += length($_)} print $_;} print "\n"; warn "\nConverted $count FASTA records in $. lines to tabular format\nTotal sequence length: $len\n\n";' enc.2593339298.genesnostop.fna  > enc.c.2593339298.tab #transform into tab-file
$ perl -e ' $col=2; while (<>) { s/\r?\n//; @F = split /\t/, $_; $len = length($F[$col])/3; print "$_\t$len\n" } warn "\nAdded column with length of column $col for $. lines.\n\n"; ' enc.c.2593339298.tab|awk -F '\t' '$4 ~ /^[0-9]+$/ { print $0 }'  |cut -f1-3 >enc.tc.2593339298.tab #select sequences which are a multiple of 3
$ perl -e ' $len=0; while(<>) { s/\r?\n//; @F=split /\t/, $_; print ">$F[0]"; if (length($F[1])) { print " $F[1]" } print "\n"; $s=$F[2]; $len+= length($s); $s=~s/.{60}(?=.)/$&\n/g; print "$s\n"; } warn "\nConverted $. tab-delimited lines to FASTA format\nTotal sequence length: $len\n\n"; ' enc.tc.2593339298.tab |sed 's/\s.*$//'> enc.tc.2593339298.genes.fna #change back to fasta and remove first space and everything after from seqids

#select sequences coding for ribosomal proteins
$ grep 'Ribosomal protein' 2593339298.cog.tab.txt|cut -f 1 >enc.2593339298.rib.list
$ perl ../extractFromFasta.pl enc.tc.2593339298.genes.fna list enc.2593339298.rib.list  > enc.tc.2593339298.ribgenes.fna

#Concatenate multiple-sequence fasta file to single long sequence
$ grep -v "^>" enc.tc.2593339298.genes.fna | awk 'BEGIN { ORS=""; print ">con_all\n" } { print }' > enc.tc.2593339298.con-all.fna
$ sed -i -e '$a\' enc.tc.2593339298.con-all.fna
$ grep -v "^>" enc.tc.2593339298.ribgenes.fna | awk 'BEGIN { ORS=""; print ">con_rib\n" } { print }' > enc.tc.2593339298.con-rib.fna
$ sed -i -e '$a\' enc.tc.2593339298.con-rib.fna

#join concatenated sequence files and a file with non-rib seqs (all sequences are doubled --> important to keep constant ratio between sequences in order to obtain correct estimations of average nucleaotide frequency which is used for downstream ENC' estimations
$ cat enc.tc.2593339298.con-rib.fna enc.tc.2593339298.con-all.fna   > ../enc.tc.2593339298.con.riball.fna
$ cat enc.tc.2593339298.con-all.fna enc.tc.2593339298.con-all.fna |sed "1s/.*/>con-rib/"  > ../enc.tc.2593339298.con.allall.fna

#remove created intermediate files
$ rm enc*

#calculate ENCprime statistics
$ cd ../
$ ENCprime-master/bin/SeqCount -c enc.tc.2593339298.con.riball.fna 2
$ ENCprime-master/bin/SeqCount -n enc.tc.2593339298.con.allall.fna 2
$ echo -en "\n\n\n0" |ENCprime-master/bin/ENCprime enc.tc.2593339298.con.riball.fna.codcnt enc.tc.2593339298.con.allall.fna.acgtfreq 11 enc.res.2593339298.txt 2
$ cat enc.res.2593339298.txt
Name Nc Ncp ScaledChi SumChi df p B_KM n_codons
con_rib: 32.6589 41.5985 0.4725 3988.8789 43 0.0000 0.5867 8442
con_all: 37.7454 47.7971 0.3040 481132.7170 43 0.0000 0.4525 1582631
Totals: 37.7269 47.7788 0.3043 484187.7488 43 0.0000 0.4526 1591073
```

The ENC prime output together with extracted ribosomal gene coding sequences were used as an input for an downstream R script, to compute &Delta;ENC' (Rocha 2004) and S (Sharp et al. 2005) and return the variables F (CUB(F)) and d.vs. (Generation.time (Vieira-Silva)).

```R
library(coRdon) #v1.1.3

# Load sequence data and the ENCprime output statistics
seq.file = list.files(pattern=glob2rx("enc.tc.*.con.riball.fna"))
ENC.file = list.files(pattern=glob2rx("enc.res.*.txt"))

i=1
seq <- readSet(file = seq.file[i])
all <- seq[2]
rib <- seq[1]
seq_mat <- codonCounts(codonTable(all))
rib_mat <- codonCounts(codonTable(rib))

P.phe.a <-sum(seq_mat[, c("TTC")])/sum(seq_mat[, c("TTC", "TTT")]) #proportion of Phe C1 codons in all genes
P.phe.r <-sum(rib_mat[, c("TTC")])/sum(rib_mat[, c("TTC", "TTT")]) #proportion of Phe C1 codons in ribosomal genes

P.ile.a <-sum(seq_mat[, c("ATC")])/sum(seq_mat[, c("ATC", "ATT")]) #proportion of Ile C1 codons in all genes
P.ile.r <-sum(rib_mat[, c("ATC")])/sum(rib_mat[, c("ATC", "ATT")]) #proportion of Ile C1 codons in ribosomal genes

P.tyr.a <-sum(seq_mat[, c("TAC")])/sum(seq_mat[, c("TAC", "TAT")]) #proportion of Tyr C1 codons in all genes
P.tyr.r <-sum(rib_mat[, c("TAC")])/sum(rib_mat[, c("TAC", "TAT")]) #proportion of Tyr C1 codons in ribosomal genes

P.asn.a <-sum(seq_mat[, c("AAC")])/sum(seq_mat[, c("AAC", "AAT")]) #proportion of Asn C1 codons in all genes
P.asn.r <-sum(rib_mat[, c("AAC")])/sum(rib_mat[, c("AAC", "AAT")]) #proportion of Asn C1 codons in ribosomal genes

# Frequency of aminoacid specific codons
freq.phe <- sum(seq_mat[, c("TTC", "TTT")])/sum(seq_mat[, c("TTC", "TTT", "ATC", "ATT", "TAC", "TAT", "AAC", "AAT")])
freq.ile <- sum(seq_mat[, c("ATC", "ATT")])/sum(seq_mat[, c("TTC", "TTT", "ATC", "ATT", "TAC", "TAT", "AAC", "AAT")])
freq.tyr <- sum(seq_mat[, c("TAC", "TAT")])/sum(seq_mat[, c("TTC", "TTT", "ATC", "ATT", "TAC", "TAT", "AAC", "AAT")])
freq.asn <- sum(seq_mat[, c("AAC", "AAT")])/sum(seq_mat[, c("TTC", "TTT", "ATC", "ATT", "TAC", "TAT", "AAC", "AAT")])

# Aminoacid specific S-values
S.phe <-log((P.phe.r*(1-P.phe.a)/P.phe.a)/(1-P.phe.r)) #aminoacid specific S-value
S.ile <-log((P.ile.r*(1-P.ile.a)/P.ile.a)/(1-P.ile.r)) #aminoacid specific S-value
S.tyr <-log((P.tyr.r*(1-P.tyr.a)/P.tyr.a)/(1-P.tyr.r)) #aminoacid specific S-value
S.asn <-log((P.asn.r*(1-P.asn.a)/P.asn.a)/(1-P.asn.r)) #aminoacid specific S-value

# S value, weighted mean of the aminoacid specific S-values
S<- S.phe*freq.phe + S.ile*freq.ile + S.tyr*freq.tyr + S.asn*freq.asn

# Delta ENCprime value
codstat <- read.table (ENC.file[i], sep=' ', header=T)
dENCp <- (codstat[2,3]-codstat[1,3])/codstat[2,3]

# Compute F and predict generation times d.vs as detailed by Vieira-Silva et al. 2010
F <- 6.747*dENCp + 1.184*S -1.438
d.vs <- (1-0.1664*(0.9726-0.7471*(1.184*S+6.747*dENCp-1.438)))^(-1/0.1664)

# GenomeID
IMG.Genome.ID<-substring(seq.file[i],8,17)

# Dataframe and display results
cub.vs=data.frame(IMG.Genome.ID,F,d.vs)
cub.vs
#IMG.Genome.ID         F     d.vs
#2593339298 0.2166112 2.389152
```

Citations:  
* Vieira-Silva, S. and Rocha, E.P.C. (2010) The Systemic Imprint of Growth and Its Uses in Ecological (Meta) Genomics. Plos Genet 6: e1000808.
* Rocha, E.P.C. (2004) Codon usage bias from tRNA’s point of view: Redundancy, specialization, and efficient decoding for translation optimization. Genome Res 14: 2279–2286.
* Sharp, P.M., Bailes, E., Grocock, R.J., Peden, J.F., and Sockett, R.E. (2005) Variation in the strength of selected codon usage bias among bacteria. Nucleic Acids Res 33: 1141–1153.
* Elek A, Kuzman M, Vlahovicek K (2020). coRdon: Codon Usage Analysis and Prediction of Gene Expressivity. R package version 1.8.0, https://github.com/BioinfoHR/coRdon. 

## Codon Usage bias (gRodon)

Growth rates predicted following a recently published model (Weissman et al 2021) were computed in R using the R package gRodon (https://github.com/jlw-ecoevo/gRodon) (see column 'Generation time (gRodon)' in Table S1). The file 2593339298.genes.fna containing sequences of all predicted genes and the file 2593339298.cog.tab.txt containing functional annotations for identifying genes encoding ribosomal proteins were used as input data.

```R
library(gRodon) #v0.0.0.9000
library(Biostrings) #v2.54.0

# Load sequence and annotation data
seq.file = list.files(pattern=glob2rx("*.genes.fna"))
ann.file = list.files(pattern=glob2rx("*.cog.tab.txt"))

# Predict generation times d.gRodon using the gRodon model
i=1
genes <- readDNAStringSet(seq.file[i])
COG <- read.table(ann.file[i], sep='\t', header=T,quote="")
rib.list <-COG[grep("*Ribosomal protein*", COG$cog_name), 1]
highly_expressed <- gsub('([0-9]+) .*', '\\1', names(genes)) %in% rib.list #edits names by removing everything after first space (including the space)
Gpred <- predictGrowth(genes, highly_expressed)
IMG.Genome.ID<-substring(seq.file[i],1,10)

# Dataframe and display results (Generation.time (gRodon)))
cub.gRodon=data.frame(IMG.Genome.ID,d.gRodon=Gpred$d)
cub.gRodon
#IMG.Genome.ID d.gRodon
#2593339298 3.360018
```

Citations:  
* Weissman, J.L., Hou, S., and Fuhrman, J.A. (2021) Estimating maximal microbial growth rates from cultures, metagenomes, and single cells via codon usage patterns. Proc Natl Acad Sci (in print)
* Elek A, Kuzman M, Vlahovicek K (2020). coRdon: Codon Usage Analysis and Prediction of Gene Expressivity. R package version 1.8.0, https://github.com/BioinfoHR/coRdon. 
* Pagès H, Aboyoun P, Gentleman R, DebRoy S (2020). Biostrings: Efficient manipulation of biological strings. R package version 2.58.0, https://bioconductor.org/packages/Biostrings.
* Henrik Bengtsson (2021). matrixStats: Functions that Apply to Rows and Columns of Matrices (and to Vectors). R package version 0.58.0. https://CRAN.R-project.org/package=matrixStats

## Species level taxonomic annotation
 
 Genome sequence data (e.g. 2593339298.fna) were copied into the directory genomes and annotated via the GTDB-Tk software (version 1.5.1, Chaumeil et al. 2019) with the following command to obtain GTDB taxonomy ([GTDB](https://gtdb.ecogenomic.org/)) annotations:

```bash

#run GTDB-Tk software (version 1.5.1)

$ gtdbtk classify_wf --genome_dir genomes --out_dir gtdb.out --cpus 80

#format gtdb-tk output

$ cd gtdb.out
$ cat gtdbtk.ar122.summary.tsv gtdbtk.bac120.summary.tsv |sort -r | sed '1d'> gtdbtk.arcbac.summary.tab
$ awk -F"\t" '{ print $1, $8,$2 }' OFS="\t" gtdbtk.arcbac.summary.tab |sed 's/classification/d__gtdb_division;p__gtdb_phylum;c__gtdb_class;o__gtdb_order;f__gtdb_family;g__gtdb_genus;s__gtdb_species/g' |sed 's/.__/\'$'\t''/g' |sed 's/;//g' |sed 's#N/A#NA#g' |cut -f1-2,4-10 >gtdb.arcbac.classification.tab
```
Genomes that based on the GTDB taxonomic annotation could not be unambiguously assigned to a single species level bin were subsequently binned at the species level via their average nucleotide identity (ANI) using reciprocal classifications via the FastANI software (version 1.32, Jain et al. 2018) followed by reformatting of the ouput file to assign species level bins with  ANI>94%.

```bash
#obtain list with genomes that could not be unambigously assigned to a single species level bin:

$ awk -F"\t" '{ print $1, $8,$2 }' OFS="\t" gtdbtk.arcbac.summary.tab |grep ';s__$' >nosp.tab #greps all genomes not assigned to the species level
$ awk -F"\t" 'NR==FNR{cnt[$3]++; next} cnt[$3]>1' nosp.tab nosp.tab |cut -f1 >../list.nosp #list with genomes not assigned at the species level 


#move genomes to be classified via FastaANI software from the directory genomes to the directory genomes.fastani and create list with path to the genome files

$ cd ../
$ LIST=$(cut -f 1 list.nosp)
$ for i in $LIST
$ do
$ mv genomes/$i.fna genomes.fastani
$ done
$ ls genomes.fastani2/* >list.fastani

#run FastANI (version 1.32)

$ fastANI --ql list.fastani --rl list.fastani -o fastani.out -t 60
$ head fastani.out
genomes.fastani/2228664027.fna    genomes.fastani/2228664027.fna    100    312    313
genomes.fastani/2228664027.fna    genomes.fastani/2716884861.fna    97.6258    189    313
genomes.fastani/2503754000.fna    genomes.fastani/2503754000.fna    100    479    479
genomes.fastani/2503754000.fna    genomes.fastani/2561511237.fna    94.8366    431    479
genomes.fastani/2503754000.fna    genomes.fastani/2739367824.fna    78.057    101    479
genomes.fastani/2508501007.fna    genomes.fastani/2508501007.fna    100    1307    1312
genomes.fastani/2508501007.fna    genomes.fastani/2546825520.fna    91.6466    1064    1312
genomes.fastani/2527291633.fna    genomes.fastani/2527291633.fna    100    319    322
genomes.fastani/2527291633.fna    genomes.fastani/2663763580.fna    81.994    105    322
genomes.fastani/2527291633.fna    genomes.fastani/2663763574.fna    81.59    173    322

#remove reciprocal matches and keep only matches ANI>94, assign species (group.pl) and remove duplicates

    #perl script group.pl to group genomes into species level bins
        #group.pl
        #!/usr/bin/perl

        use strict;

        my $gc = 1; # group counter
            my %species; # hash containing group numbers for each element
                my @pairs;  # array of arrays containing pairs

        while(<>) {
        chomp;
        my ($a,$b) = split /\t/;

        $species{$a} = $gc++ unless (defined($species{$a}));

        $species{$b} = $species{$a};
        push @{ $pairs[$species{$a}] }, [ $a, $b ];
        };

        END {
        for my $g (keys @pairs) {
        for my $p (@{ $pairs[$g] }) {
        printf "%s\t%s\tfastani_species%02i\n", @$p[0], @$p[1], $g;
        }
        };
        }

$ awk '$1!=$2' fastani.out |awk '$3 > 94' |cut -f1,2 | sed 's#genomes.fastani/##g' |sed 's/.fna//g' |perl group.pl |cut -f1,3 |sort|uniq >fastani.bins
$ head fastani.bins
2228664027    fastani_species01
2503754000    fastani_species02
2561511237    fastani_species02
2596583567    fastani_species03
2608642245    fastani_species04
2615840672    fastani_species05
2651869507    fastani_species03
2651869509    fastani_species06
2651869521    fastani_species03
2651869522    fastani_species06
```
Citations:  
* Parks, D.H., et al. (2022)." GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy." Nucleic Acids Research
* Chaumeil, P.-A, et al. (2019). "GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database." Bioinformatics, btz848: https://doi.org/10.1093/bioinformatics/btz848. 
* Jain, C., Rodriguez-R, L. M., Phillippy, A. M., Konstantinidis, K. T., and Aluru, S. (2018). High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. Nat. Commun. 9, 5114. doi:10.1038/s41467-018-07641-9.
