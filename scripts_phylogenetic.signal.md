# Scripts and comand lines for the computation of phylogenetic signals and mantel correlograms

### by Sara Beier with support of Johannes Werner

This file contains R scripts to compute phylogenetic signals and matel correlograms of genomic traits from the genomic traits values provided in Table S1 (Beier et al, doi: https://doi.org/10.3389/fmicb.2022.985216), with the only exception of the 16s rRNA gene copy number (RRN). Due to detected biases of RRN values available via the [JGI/IMG platform](https://img.jgi.doe.gov/) (see Figure S3) we use in this case trait values and phylogeny of the the curated Ribosomal RNA Operon Copy Number Database ([rrnDB](https://rrndb.umms.med.umich.edu/)).


## rrnDB phylogeny

16s rRNA sequences (rrnDB-5.7_16S_rRNA.fasta) from entries stored in the [rrnDB](https://rrndb.umms.med.umich.edu/) (Stoddard et al. 2015) were alligned after removing duplicate sequences using the MUSCLE software (version 3.8.1551, Edgar 2004). A phylogenetic tree was computed using the FastTree software  (version 2.1.10, Price et al. 2010).

```bash
# Change fasta to tab
perl -e ' $count=0; $len=0; while(<>) { s/\r?\n//; s/\t/ /g; if (s/^>//) { if ($. != 1) { print "\n" } s/ |$/\t/; $count++; $_ .= "\t"; } else { s/ //g; $len += length($_) } print $_; } print "\n"; warn "\nConverted $count FASTA records in $. lines to tabular format\nTotal sequence length: $len\n\n"; ' rrnDB-5.7_16S_rRNA.fasta >rrnDB-5.7_16S_rRNA.tab

# Remove duplicate sequences
perl -e '$column = 2; $unique=0; while(<>) {s/\r?\n//; @F=split /\t/, $_; if (! ($save{$F[$column]}++)) {print "$_\n"; $unique++}} warn "\nChose $unique unique lines out of $. total lines.\nRemoved duplicates in column $column.\n\n"' rrnDB-5.7_16S_rRNA.tab |cut -f 2,5 -d '|' | sed 's/|/\'$'\t''/g' |cut -f1,3 > unique.tab

# Keep only one of duplicates in column 1, add one empty column (necessary to retransform into fasta)
awk '!a[$1]++' unique.tab  |sed 's/\'$'\t''/\'$'\t''\'$'\t''/g' > rrnDB-5.7_16S_rRNA.dedup.tab

# Change tab to fasta
perl -e ' $len=0; while(<>) { s/\r?\n//; @F=split /\t/, $_; print ">$F[0]"; if (length($F[1])) { print " $F[1]" } print "\n"; $s=$F[2]; $len+= length($s); $s=~s/.{60}(?=.)/$&\n/g; print "$s\n"; } warn "\nConverted $. tab-delimited lines to FASTA format\nTotal sequence length: $len\n\n"; ' rrnDB-5.7_16S_rRNA.dedup.tab > rrnDB-5.7_16S_rRNA.dedup.fasta

# Sequence alignement using the muscle aligner
muscle -in rrnDB-5.7_16S_rRNA.dedup.fasta -out muscle.rrnDB.fasta -maxiters 2

# Compute phylogenetic tree with the FastTree software
FastTree -gamma -gtr -nt -log rrndb_fasttree.log muscle.rrnDB.fasta

# Format output tree
cat rrndb_fasttree.tree | sed '1,/Total time:/d' > rrndb_fasttree.format.tree #remove information lines from fasttree output tree (otherwise tre cannot be read by phylools (R, physig2021_rrndb_fasttree.R))
```
Citations:  
* Stoddard, S. F., Smith, B. J., Hein, R., Roller, B. R. K., and Schmidt, T. M. (2015). rrnDB: improved tools for interpreting rRNA gene abundance in bacteria and archaea and a new foundation for future development. Nucleic Acids Res. 43, D593–D598. doi:10.1093/nar/gku1201..
* Edgar, R. C. (2004). MUSCLE: a multiple sequence alignment method with reduced time and space complexity. BMC Bioinformatics 5, 1–19. doi:10.1186/1471-2105-5-113.
* Price, M. N., Dehal, P. S., and Arkin, A. P. (2010). FastTree 2-Approximately Maximum-Likelihood Trees for Large Alignments. Plos One 5. doi:10.1371/journal.pone.0009490.

## Phylogenetic signals

The phylogenetic signal of each trait was assessed via Pagel's lambda (Pagel 1999) and Blomberg's K (Blomberg et al. 2003) statistics.  
The trait 'Generation time' was for all downstream statistics  log(x) transformed and the trait %HGT was log(x+0.01) transformed. Visual inspection of the data  revealed a better fit to normal distribution after these transformation steps. It was not possible to apply transformations for the highly skewed distributions of the number of16s rRNA gene copies (RRN_rrnDB) or the number of prophages that would have resulted in a better fit to normal distribution. 
For all analyses except for the trait RRN, the phylogeny of the prokaryotic PICRUSt2 default phylogenetic tree (Douglas et a. 2020: pro_ref.tre) was used to infer phylogenetic signals of genomic traits among the PICRUSt2 reference genomes. To infer phylogenetic signals for RRN we used the trait table available viathe [rrnDB](https://rrndb.umms.med.umich.edu/)  and a phylogenetic tree that was computed from the corresponding 16s rRNA gene sequence data.

```R
library(phytools) #v0.7.20

# Load data
pic.tre <-read.newick("pro_ref.tre") #Prokaryote reference tree from the PICRUSt2 softaware
rrn.tre <- read.newick ("rrndb_fasttree.format.tree") # rrnDB tree computed as detailed above
gtraits <- read.table ("TableS1.tsv", header=T, sep='\t', fill=T)[,c(2,18,20:21,23:24, 26:30)]
gtraits <-gtraits[gtraits$prophages<26 & gtraits$X.HGT <= 65,] #remove outlayers: prophages<26 and %HGT<=65
rrnDBtraits <- read.csv ("rrnDB-5.7.tsv", sep='\t', header=T, fill=T)[-c(1:215),c(1,12)] #rrnDB trait table, exclude first rows, no sequence data available

# Select trait for analyses 
trait <-gtraits[,2] # Genome size
#trait <-gtraits[,4] # %GC
#trait <-log(gtraits[,5]+0.01) # %HGT (log(x+0.01) transformation!)
#trait <-gtraits[,6] # CUB (F)
#trait <-log(gtraits[,7]) # Generation time (gRodon) (log(x) transformation!)
#trait <-gtraits[,8] # Gene duplication
#trait <-gtraits[,9] # Gene richness
#trait <-gtraits[,10] # %TF
#trait <-gtraits[,11] # Prophages
names(trait) <-gtraits[,1]

#trait <-rrnDBtraits[,2] # RRN
#names(trait) <-rrnDBtraits[,1]

# Test for phylogenetic signal using Blomberg's K statistics
Kstats <-phylosig(pic.tre, na.omit(trait), method="K", test=TRUE, nsim=1001) #for traits from Table S1
#Kstats <-phylosig(rrn.tre, na.omit(trait), method="K", test=TRUE, nsim=1001) #for RRN

# Test for phylogenetic signal using Pagel's lamda statistics
Lstats <-phylosig(pic.tre, na.omit(trait), method="lambda", test=TRUE, nsim=1001) #for traits from Table S1
#Lstats <-phylosig(rrn.tre, na.omit(trait), method="lambda", test=TRUE, nsim=1001) #for RRN
```
Citations:  
* Pagel, M. (1999) Inferring the historical patterns of biological evolution. Nature 401: 877–884.
* Blomberg, S.P., Garland, T., and Ives, A.R. (2003) Testing for phylogenetic signal in comparative data: Behavioral traits are more labile. Evolution 57: 717–745.
* Douglas, G.M., Maffei, V.J., Zaneveld, J.R., Yurgel, S.N., Brown, J.R., Taylor, C.M., et al. (2020) PICRUSt2 for prediction of metagenome functions. Nat Biotechnol 38: 685–688.


## Mantel correlograms

We computed phylogenetic Mantel correlograms similar as detailed elsewhere (Diniz-Filho et al. 2010, Dini-Andreote et al. 2015) to test for significant positive correlations between phylogenetic distances and trait distances at the following phylogenetic distance classes: 0-0.25 / 0.25-0.5 / 0.5-0.75 / 0.75-1 / 1-1.5 / 1.5-2 / 2-2.5/ 2.5-3. For the calculation of Mantel correlograms, the dataset was reduced to 10000 randomly selected genomes in order to reduce computation time and memory demand. The Mantel correlations were tested via the Spearman rank-order correlation as it was not possible to obtain normal distributed trait values in all cases. 
For all analyses except for the trait RRN, the phylogeny of the prokaryotic PICRUSt2 default phylogenetic tree (Douglas et a. 2020: pro_ref.tre) was used to infer phylogenetic signals of genomic traits among the PICRUSt2 reference genomes. To infer phylogenetic signals for RRN we used the trait table available via [rrnDB](https://rrndb.umms.med.umich.edu/) and a phylogenetic tree that was computed from the corresponding 16s rRNA gene sequence data.

```R
library(phytools) #v0.7.20
library(mpmcorrelogram) #v0.1.4
'%!in%' <- function(x,y)!('%in%'(x,y)) #define 'not in' function

# Load data
pic.tre <-read.newick("pro_ref.tre") #Prokaryote reference tree from the PICRUSt2 softaware
rrn.tre <- read.newick ("rrndb_fasttree.format.tree") # rrnDB tree computed as detailed above
gtraits <- read.table ("TableS1.tsv", header=T, sep='\t', fill=T)[,c(2,18,20:21,23:24, 26:30)]
gtraits <-gtraits[gtraits$prophages<26 & gtraits$X.HGT <= 65,] #remove outlayers: prophages<26 and %HGT<=65

rrnDBtraits <- read.csv ("rrnDB-5.7.tsv", sep='\t', header=T, fill=T)[-c(1:215),c(1,12)] #rrnDB trait table, exclude first rows, no sequence data available

# Set parameters
n = 10000 # Set n for subsampling
perm = 200 # Set number of permutations for Mantel correlograms
breaks=c(0,0.25,0.5,0.75,1,1.5,2,2.5,3) # Set breaks for Mantel correlograms

# Select trait for analyses 
trait <-gtraits[,2] # Genome size
#trait <-gtraits[,4] # %GC
#trait <-log(gtraits[,5]+0.01) # %HGT (log(x+0.01) transformation!)
#trait <-gtraits[,6] # CUB (F)
#trait <-log(gtraits[,7]) # Generation time (gRodon) (log(x) transformation!)
#trait <-gtraits[,8] # Gene duplication
#trait <-gtraits[,9] # Gene richness
#trait <-gtraits[,10] # %TF
#trait <-gtraits[,11] # Prophages
names(trait) <-gtraits[,1]

#trait <-rrnDBtraits[,2] #for RRN
#names(trait) <-rrnDBtraits[,1]

# Randomly select 10000 entries and prune trees (traits from Table 1)
set.seed(2022)
gtraits <- na.omit(gtraits) # Remove lines with NA from trait table
tips<-sample(pic.tre$tip.label) #all species in tree
tips.n <- sample(gtraits$PICRUST.ID,n) #randomly select 1000- species
tips.notn <- tips[tips %!in% tips.n] #not selected species
tre.sub <-fancyTree(pic.tre,type="droptip",tip=tips.notn) #prune tree by removing not selected species

# Randomly select 10000 entries and prune tree (RRN)
set.seed(2022)
tips<-sample(rrn.tre$tip.label) #all species in tree
tips.n <- sample(tips,n) #randomly select 1000- species
tips.notn <- tips[tips %!in% tips.n] #not selected species
tre.rrn.sub <-fancyTree(rrn.tre,type="droptip",tip=tips.notn) #prune tree by removing not selected species


# Compute pairwise phylogenetic distances (traits from Table 1)
dist.p <- cophenetic(tre.sub) 
#dist.p <- cophenetic(tre.rrn.sub) #RRN

# Compute pairwise distances for trait values
dist.t <- dist(trait[names(trait) %in% tre.sub$tip.label], method = 'manhattan')
#dist.t <- dist(trait[names(trait) %in% tre.sub.rrn$tip.label], method = 'manhattan') #RRN

# Mantel correlogram (traits from Table 1)
set.seed(2022)
man <- mpmcorrelogram(dist.t, dist.p, method="spearman", permutations=perm, breaks=breaks,plot = FALSE, print = FALSE) # Mantel correlogram
res <- do.call(cbind.data.frame, man[c(6,2,4,5)])
res$pd <- breaks[2:length(breaks)]
res$trait <-colnames(gtraits[trait])
res
write.table(res, "mpm.2.tab",sep = '\t', row.names = FALSE) #save as mpm.2.tab, mpm.4.tab,... for gtraits[,2], gtraits[,4],...
```

Citations:  
* Dini-Andreote, F., Stegen, J.C., Elsas, J.D. van, and Salles, J.F. (2015) Disentangling mechanisms that mediate the balance between stochastic and deterministic processes in microbial succession. Proc Natl Acad Sci 112: E1326–E1332.
* Diniz-Filho, J.A.F., Terribile, L.C., Cruz, M.J.R. da, and Vieira, L.C.G. (2010) Hidden patterns of phylogenetic non-stationarity overwhelm comparative analyses of niche conservatism and divergence. Glob Ecol Biogeogr 19: 916–926.
* Douglas, G.M., Maffei, V.J., Zaneveld, J.R., Yurgel, S.N., Brown, J.R., Taylor, C.M., et al. (2020) PICRUSt2 for prediction of metagenome functions. Nat Biotechnol 38: 685–688.
