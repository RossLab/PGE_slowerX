# 6. Effective population size

To test our prediction that the effective population sizes of the X and autosomes should be the same under PGE, i.e. NeX/NeA = 1, I looked at neutral diversity. As a measure of neutral diversity, I estimated heterozygosity at fourfold degenerate sites using ANGSD.

To generate degeneracy annotations I used the script partitioncds.py from Mackintosh et al. (2022) - see 04_divergence.md. One of the intermediate files I produce is a *.4x.sites.bed file, which looks like this:

```
head Ling.augustus.genes.4x.sites.fixed.bed
scaffold_0_LN_i_129109_RC_i_216858_XC_f_1.000000	12519	12520	g21688.t1	1	AUGUSTUS	mRNA	.	g21688.t1
scaffold_0_LN_i_129109_RC_i_216858_XC_f_1.000000	12510	12511	g21688.t1	1	AUGUSTUS	mRNA	.	g21688.t1
scaffold_0_LN_i_129109_RC_i_216858_XC_f_1.000000	12501	12502	g21688.t1	1	AUGUSTUS	mRNA	.	g21688.t1
scaffold_0_LN_i_129109_RC_i_216858_XC_f_1.000000	12495	12496	g21688.t1	1	AUGUSTUS	mRNA	.	g21688.t1
```

And can be subsetted to get a list of sites that can be passed to ANGSD:

```
cut -f1,3 Ling.augustus.genes.4x.sites.fixed.bed > Ling.XX.4x.sites
head Ling.XX.4x.sites
scaffold_0_LN_i_129109_RC_i_216858_XC_f_1.000000	12520
scaffold_0_LN_i_129109_RC_i_216858_XC_f_1.000000	12511
scaffold_0_LN_i_129109_RC_i_216858_XC_f_1.000000	12502
scaffold_0_LN_i_129109_RC_i_216858_XC_f_1.000000	12496
```

ANGSD also needs a list of bam files, i.e.

```
Ling.XX.bam.filelist
L19_EKDN230015189-1A_HYCGGDSX5_L1.Lycoriella_ingenua_scaffolds.rascaf.i1.fixed.fa.pic.rmdup.rg.sorted.bam
L20_EKDN230015190-1A_HYCGGDSX5_L1.Lycoriella_ingenua_scaffolds.rascaf.i1.fixed.fa.pic.rmdup.rg.sorted.bam
L21_combined.Lycoriella_ingenua_scaffolds.rascaf.i1.fixed.fa.pic.rmdup.rg.sorted.bam
L22_EKDN230015192-1A_HYCGGDSX5_L1.Lycoriella_ingenua_scaffolds.rascaf.i1.fixed.fa.pic.rmdup.rg.sorted.bam
L23_EKDN230015193-1A_HYCGGDSX5_L1.Lycoriella_ingenua_scaffolds.rascaf.i1.fixed.fa.pic.rmdup.rg.sorted.bam
```

ANGSD will then take these ingredients and generate a site frequency spectrum, from which it can be used to calculate theta (heterozygosity):

```
angsd sites index Ling.XX.4x.sites
angsd -bam Ling.XX.bam.filelist -sites Ling.XX.4x.sites -doSaf 1 -anc $GENOME -GL 1 -P 24 -out Ling.XX.angsd 
#for folded
realSFS Ling.XX.angsd.saf.idx -P 8 -fold 1 > Ling.XX.angsd.sfs
realSFS saf2theta Ling.XX.angsd.saf.idx -outname Ling.XX.angsd -sfs Ling.XX.angsd.sfs -fold 1
#Estimate for every Chromosome/scaffold
thetaStat do_stat Ling.XX.angsd.thetas.idx
```

I then calculated NeX/NeA in RStudio:

```
### Lycoriela ingenua

# angsd output
lingNe <- read.table('Ling.XX.angsd.thetas.idx.pestPG.fixed', header=F, stringsAsFactors=F)
lingNe <- lingNe[,c(2,4,14)]
colnames(lingNe) <- c('chr', 'tW', 'nSites')
# sex linkage
LingSex <- read.table('../sex_linkage/assignments/ling_sexlinkage.fixed.txt', header=F, stringsAsFactors=F)
colnames(LingSex) <- c('chr', 'linkage')

lingNe <- merge(lingNe, LingSex, by=c('chr'))

lingNeA <- lingNe[which(lingNe$linkage=="A"),]
lingNeX <- lingNe[which(lingNe$linkage=="X"),]

lingNeA <- lingNeA[which(lingNeA$tW > 0),]
lingNeX <- lingNeX[which(lingNeX$tW > 0),]

# because each scaffold doesn't have the same number of sites, I need to 'weigh' theta for each scaffold by the number of sites
sum(lingNeA$tW) # 8349.751
sum(lingNeA$nSites) # 3938285
sum(lingNeX$tW) # 3441.011
sum(lingNeX$nSites) # 1659603
8349.751/3938285 # 0.002120149
3441.011/1659603 # 0.002073394
0.002073394/0.002120149 # 0.9779473
```




















