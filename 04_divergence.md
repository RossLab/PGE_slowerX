# 4. Divergence

### 4.1. Prep - identifying species pairs

Before I go ahead and calculate _dN/dS_, I wanted to (i) identify which species pairs to use, (ii) check divergence (Dxy) between species within each pair, and (iii) identify orthologous genes within pairs from which to calculate _dN/dS_.

To identify which species pairs to use, I generated a phylogeny using the BUSCO output from when I generated the genome assemblies (I had to run BUSCO on all the other genomes too, using the insecta database). Dilophus (Bibionidae) will be the outgroup. I used BUSCO phylogenomics:

```
python ~/software/BUSCO_phylogenomics/BUSCO_phylogenomics.py -i BUSCO_insecta -o BUSCO_phylo_insecta_80_iqtree -t 8 -psc 80.0 --gene_tree_program iqtree
```

And then IQTree on the supermatrix, using the parameters Q.insect+F+I+R5:

```
iqtree -s SUPERMATRIX.fasta -alrt 1000 -bb 1000 -nt AUTO -ntmax 8 -m Q.insect+F+I+R5 
```

I loaded the phylogeny into FigTree and assigned the following species pairs. In all cases, I decided that I would use the species with the most contiguous genome assembly as the focal/ingroup, and would map reads from the other species to calulcate dN/dS:

| Focal species | Outgroup species |
| - | - |
| Phytosciara flavipes | Bradysia pectoralis |
| Bradysia coprophila | Bradysia odoriphaga |
| Bradysia confinis | Bradysia desolata |
| Lycoriella ingenua | Lycoriella agraria |
| Contarinia nasturtii | Contarinia rumicis |
| Sitodiplosis mosellana | Aphidoletes aphidimyza |
| Mayetiola destructor | Mayetiola hordei |
| Dilophus febrilis | Dilophus femoratus |

Next, I ran alignments between the genomes with minimap2:

```
minimap2 -t 30 -c --secondary=no -ax asm20 $FOCAL*.fa $OUTGROUP*.fa | samtools view -bo $FOCAL.$OUTGROUP.bam && samtools sort -o $FOCAL.$OUTGROUP.sorted.bam $FOCAL.$OUTGROUP.bam && rm $FOCAL.$OUTGROUP.bam
bedtools bamtobed -ed -i $FOCAL.$OUTGROUP.sorted.bam > $FOCAL.$OUTGROUP.sorted.bed && rm $FOCAL.$OUTGROUP.sorted.bam
```

And calculated Dxy as the 1 - alignment length / number of mismatches:

```
aln <- read.table('dxy/Smos.Aaph.sorted.bed', header=F, stringsAsFactors=F)
colnames(aln) <- c('ctg_ref', 'aln_start', 'aln_end', 'ctg_query', 'N_mismatches', 'orientation')
aln$aln_len <- aln$aln_end-aln$aln_start
# filter out super short alignments
nrow(aln)
aln <- aln[which(aln$aln_len > 1000),]
nrow(aln)
# sum aln length and NM
# Dxy = 1-PI, i.e. NM/aln
dxy <- sum(aln$N_mismatches)/sum(aln$aln_len)
dxy
```

Then, I ran Orthofinder to identify single-copy orthologs between species within each pair, which I would use to calculate _dN/dS_. I excluded genes that were sex-linked in one species but autosomal in the other from further analysis. Where sets of predicted genes had multiple isoforms for each gene (this was the case for all BRAKER predictions), I took the longest isoform for each gene.

```
OUT=Smos_v_Aaph
mv Smos_genes.fasta genes/
mv Aaph_genes.fasta genes/
orthofinder -o $OUT -n $OUT -t 12 -d -f genes/
mv $OUT/Results_$OUT/Orthogroups/Orthogroups.tsv $OUT.OGs.tsv
mv $OUT/Results_$OUT/Orthogroups/Orthogroups_SingleCopyOrthologues.txt $OUT.SCOs.tsv
```

### 4.2. Site counts: synonymous and nonsynonymous sites per gene

In order to calculate _dN/dS_, I need to run degeneracy annotations for each species. This is to get the number of synonymous and non-synonymous _sites_ per gene, i.e. the denominators in _dN/dS_. Remember that _dN/dS_ is the ratio of the number of nonsynonymous substitutions per non-synonymous site to the number of synonymous substitutions per synonymous site.

I will use Dom's script for Alex's paper to annotate degeneracy: https://github.com/A-J-F-Mackintosh/Mackintosh_et_al_2022_Ipod/tree/main

It requires three inputs:
1. A genome (note it won't work if there are whitespaces in fa headers)
2. An annotation bed file
3. a VCF

Note that this script takes a long time to run on highly fragmented genomes. The main limiting factor is how long it takes to parse the VCF, which is proportional to the number of contigs, although the number of variants in the VCF also affects the speed so it could be worth applying some harder filters if time matters. I've noticed that a genome with 4 contigs takes minutes, while one with >100k contigs can take over two weeks.

It is not strictly necessary to use a VCF to annotate codon degeneracy, but unfortunately this script only works if you supply a VCF that contains variants (an empty VCF won't work). The purpose of the VCF is to provide a more accurate degeneracy annotation. This is because for any given codon, whether there is a SNP or a het site matters for its degeneracy: for example, imagine a codon where the third position is a fourfold degenerate site. If there's a SNP at the first site then it may no longer be fourfold degenerate. Using a vcf is unlikely to change things much - I tested this by suppling a full VCF and then only the first ~20k variants (or something like that...) and it only affected a minority of genes, and the difference was subtle.

For a couple of the species I'm analysing I have population data, but for most I do not, so I'll simply align an Illumina library for each focal species to their respective genomes (or PacBio if no Illumina is available) and call variants to generate the VCF. 

Align with bwa mem:

```
GENOME=*.fa
bwa index ${base}.*.fa
bwa mem -t 32 ${base}.*.fa ${base}_1.trimmed.fastq.gz ${base}_2.trimmed.fastq.gz | \
samtools view -b - > ${base}.bwa.bam && \
samtools sort -@ 8 -o ${base}.bwa.sorted.bam ${base}.bwa.bam
```

Process alignments with picard tools (sort, mark and remove duplicates, add read groups):
```
# sort with picard
picard SortSam I=${base}.bwa.sorted.bam SORT_ORDER=coordinate O=${base}.bwa.pic.sorted.bam && \
# mark and remove duplicates
samtools index ${base}.bwa.pic.sorted.bam && \
java -Xmx20g -jar /ceph/users/rbaird/software/picard/build/libs/picard.jar \
MarkDuplicates MAX_RECORDS_IN_RAM=400000 INPUT=${base}.bwa.pic.sorted.bam OUTPUT=${base}.bwa.pic.rmdup.sorted.bam M=${base}.met.txt REMOVE_DUPLICATES=true --READ_NAME_REGEX && \
# assign unique read groups
picard AddOrReplaceReadGroups INPUT=${base}.bwa.pic.rmdup.sorted.bam
OUTPUT=${base}.bwa.pic.rmdup.rg.sorted.bam RGID=${base} RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${base}
```

Call, genotype and filter variants with gatk-4
```
/ceph/users/rbaird/software/gatk/gatk CreateSequenceDictionary -R $GENOME
# call variants
samtools index ${base}.bwa.pic.rmdup.rg.sorted.bam && \
/ceph/users/rbaird/software/gatk/gatk --java-options "-Xmx4g" HaplotypeCaller \
-R $GENOME \
-I ${base}.bwa.pic.rmdup.rg.sorted.bam \
-O ${base}.gatk.vcf.gz \
-ERC GVCF && \
# genotype variants
/ceph/users/rbaird/software/gatk/gatk --java-options "-Xmx4g" GenotypeGVCFs \
-R $GENOME \
-V ${base}.gatk.vcf.gz \
-O ${base}.gatk.g.vcf.gz
/ceph/users/rbaird/software/gatk/gatk SelectVariants \
-R $GENOME \
-V ${base}.gatk.g.vcf.gz \
-O ${base}.gatk.g.filt.vcf.gz \
-select "QD > 2.0" \
-select "FS < 60.0" \
-select "MQ > 40.0"
```

The bed file that partitioncds.py requires looks like this:
```
CM059996.1	618	891	g11661.t1	1	-
CM059996.1	1413	1548	g11661.t1	0.96	-
CM059996.1	3671	3842	g11661.t1	0.71	-
CM059996.1	6509	6558	g11662.t1	1	+
CM059996.1	6642	6767	g11662.t1	1	+
```

Where columns 2 and 3 are CDS coordinates. This can be generated from a gff3/gtf using GenomeTools (and a bit of magic):
```
# this works with an augustus.hints.gff3 output from braker2, it may need to be adjusted if using an annotation generated by other means... Look into the tool AGAT for fixing GFFs.
for file in $(ls *.gff3)
do
	base=$(basename $file ".gff3")
	cat ${base}.gff3 | gt gff3 -sort -tidy -retainids -fixregionboundaries > ${base}.tidy.gff3
	gff2bed < ${base}.tidy.gff3 > ${base}.tidy.gt.bed
	cat ${base}.tidy.gt.bed | awk '$8=="CDS"' OFS="\t" | cut -f 1,2,3,4,5,6 | sort -k1,1V -k2,2n | awk 'BEGIN{FS=OFS="\t"} {sub(/.CDS.*/,"",$4)} 1' > ${base}.tidy.gt.CDS.bed
done
```

Now generate the degeneracy annotation:
```
GENOME=*.fa
BEDFILE=*.gt.CDS.bed
VCF=*.filt.vcf.gz
ANNO=*.anno.bed

~/software/Mackintosh_et_al_2022_Ipod/partitioncds.py -f $GENOME -v $VCF -b $BEDFILE
```

The output gimble.cds.bed contains the degeneracy annotation in the following format:
```
CM059996.1	3840	3841	M_0	2	-
CM059996.1	3839	3840	M_0	3	-
CM059996.1	3838	3839	F_0	1	-
CM059996.1	3837	3838	F_0	2	-
CM059996.1	3836	3837	F_2	3	-
```

Where the fifth column contains the degeneracy annotation - e.g. M_0 means Methionine_zero-fold. I subsetted the 0x, 2x, 3x and 4x degenerate sites:

```
### subset 0x, 2x, 3x and 4x degenerate sites:
grep "_0" gimble.cds.bed > gimble.cds.0xsites.bed 
grep "_2" gimble.cds.bed > gimble.cds.2xsites.bed 
grep "_3" gimble.cds.bed > gimble.cds.3xsites.bed 
grep "_4" gimble.cds.bed > gimble.cds.4xsites.bed 
gzip gimble.cds.bed
```

...and then intersected this with a bed file to include a gene column (the script doesn't output that unfortunately):

```
gff2bed < anno.gff3 > anno.bed
ANNO=anno.bed
### intersect degeneracy annotation with bedfile
bedtools intersect -a *.anno.bed -b gimble.cds.0xsites.bed > genes.0x.sites.bed && rm gimble.cds.0xsites.bed
grep "mRNA" genes.0x.sites.bed | sed 's/ID=//g' | sed 's/;.*//g' | cut -f1,2,3,10 > genes.0x.sites.fixed.bed && gzip genes.0x.sites.fixed.bed && rm genes.0x.sites.bed
bedtools intersect -a *.anno.bed -b gimble.cds.2xsites.bed > genes.2x.sites.bed && rm gimble.cds.2xsites.bed
grep "mRNA" genes.2x.sites.bed | sed 's/ID=//g' | sed 's/;.*//g' | cut -f1,2,3,10 > genes.2x.sites.fixed.bed && gzip genes.2x.sites.fixed.bed && rm genes.2x.sites.bed
bedtools intersect -a *.anno.bed -b gimble.cds.3xsites.bed > genes.3x.sites.bed && rm gimble.cds.3xsites.bed
grep "mRNA" genes.3x.sites.bed | sed 's/ID=//g' | sed 's/;.*//g' | cut -f1,2,3,10 > genes.3x.sites.fixed.bed && gzip genes.3x.sites.fixed.bed && rm genes.3x.sites.bed
bedtools intersect -a *.anno.bed -b gimble.cds.4xsites.bed > genes.4x.sites.bed && rm gimble.cds.4xsites.bed
grep "mRNA" genes.4x.sites.bed | sed 's/ID=//g' | sed 's/;.*//g' | cut -f1,2,3,10 > genes.4x.sites.fixed.bed && gzip genes.4x.sites.fixed.bed && rm genes.4x.sites.bed
# get degen site counts
Rscript  ~/scripts/R_scripts/get_degen_counts.R
```

The R script will count the number of synonymous and non-synonymous sites per gene (a two-fold site is counted as 0.5 syn and 0.5 nonsyn, a 3x site as 0.66 syn and 0.22 nonsyn):
```
# R

library(data.table)
library(dplyr)

degen_0_sites <- read.table('genes.0x.sites.fixed.bed.gz', header=F, stringsAsFactors=F)
colnames(degen_0_sites) <- c('scaf', 'bed_start', 'bed_end', 'transcript')
degen_2_sites <- read.table('genes.2x.sites.fixed.bed.gz', header=F, stringsAsFactors=F)
colnames(degen_2_sites) <- c('scaf', 'bed_start', 'bed_end', 'transcript')
degen_3_sites <- read.table('genes.3x.sites.fixed.bed.gz', header=F, stringsAsFactors=F)
colnames(degen_3_sites) <- c('scaf', 'bed_start', 'bed_end', 'transcript')
degen_4_sites <- read.table('genes.4x.sites.fixed.bed.gz', header=F, stringsAsFactors=F)
colnames(degen_4_sites) <- c('scaf', 'bed_start', 'bed_end', 'transcript')

degen_0_counts <- degen_0_sites %>% count(transcript)
colnames(degen_0_counts) <- c('transcript', '0x_sites')
degen_2_counts <- degen_2_sites %>% count(transcript)
colnames(degen_2_counts) <- c('transcript', '2x_sites')
degen_3_counts <- degen_3_sites %>% count(transcript)
colnames(degen_3_counts) <- c('transcript', '3x_sites')
degen_4_counts <- degen_4_sites %>% count(transcript)
colnames(degen_4_counts) <- c('transcript', '4x_sites')

degen_counts_02 <- merge(degen_0_counts, degen_2_counts, by=c('transcript'), all=T)
degen_counts_023 <- merge(degen_counts_02, degen_3_counts, by=c('transcript'), all=T)
all_degen_counts <- merge(degen_counts_023, degen_4_counts, by=c('transcript'), all=T)

# get N syn vs nonsyn sites
all_degen_counts$N_non <- (all_degen_counts$`0x_sites`)+(all_degen_counts$`2x_sites`*0.5)+(all_degen_counts$`3x_sites`*0.33)
all_degen_counts$N_syn <- (all_degen_counts$`4x_sites`)+(all_degen_counts$`2x_sites`*0.5)+(all_degen_counts$`3x_sites`*0.66)
nrow(all_degen_counts)
# write table
write.table(all_degen_counts, file='degen_site_counts.tsv', row.names=F, col.names=T, quote=F, sep='\t')
```

### 4.3. Variant counts: alignments and variant calling

Now comes the numerators: the nonsynonymous and synonymous substitutions. 

I generally use bwa mem to align Illumina data from an outgroup to each focal group where possible. I find bwa gives better mapping rates than Bowtie2, although it's a little slower. If using PacBio (HiFi), I use minimap2. See above for how I decided which outgroups to use for each focal group.

E.g.:
```
FOCAL=Mdes
OUTGR=Mhor
GENOME=Mdes.ncbi.msk
R1=Mhor_1.trimmed.fastq.gz
R2=Mhor_2.trimmed.fastq.gz
bwa index $GENOME
bwa mem -t 32 $GENOME $R1 $R2 | samtools view -b - > $FOCAL.$OUTGR.bwa.bam && samtools sort -@ 8 -o $FOCAL.$OUTGR.bwa.sorted.bam $FOCAL.$OUTGR.bwa.bam && rm $FOCAL.$OUTGR.bwa.bam && samtools index $FOCAL.$OUTGR.bwa.sorted.bam
# Check mapping rate
samtools flagstats $FOCAL.$OUTGR.bwa.sorted.bam > $FOCAL.$OUTGR.bwa.sorted.bam.flagstats
```

Then process bam files with Picardtools - sort, mark/remove duplicates, add read groups. And then call, genotype and filter variants with gatk-4. See above - I used the same commands as before.

Finally, variant files need to be annotated for their effects. I use SnpEff/SnpSift for this.

The snpEff.config file needs a line added to it so that you can call a species, e.g.:
```
nano SnpEff/snpEff.config
## Scroll down and add the following:
# Mayetiola destructor genome (GCA_000149185.1)
Mayetiola_destructor.genome : Mayetiola_destructor
```

Then the genome and annotation needs to go in SnpEff/data/ before a SnpEff database can be built. I can only get this to work with gtfs, and only if I convert my gff3s to gtfs (it doesn't like the gtf that braker2 outputs). So I convert to gtf with gffread first.
```
SPECIES=Mayetiola_destructor
SP_PREFIX=Mdes
gffread $SP_PREFIX*.gff3 -T -o $SP_PREFIX.gtf
mkdir -p ~/.conda/envs/SnpEff/SnpEff/data/$SPECIES
cp $SP_PREFIX*.fa ~/.conda/envs/SnpEff/SnpEff/data/$SPECIES
cp $SP_PREFIX.gtf ~/.conda/envs/SnpEff/SnpEff/data/$SPECIES
mv ~/.conda/envs/SnpEff/SnpEff/data/$SPECIES/*.fa ~/.conda/envs/SnpEff/SnpEff/data/$SPECIES/sequences.fa
mv ~/.conda/envs/SnpEff/SnpEff/data/$SPECIES/*.gtf ~/.conda/envs/SnpEff/SnpEff/data/$SPECIES/genes.gtf
# Build database (check if output snpEffectPredictor.bin was produced)
java -Xmx20g -jar ~/.conda/envs/SnpEff/SnpEff/snpEff.jar build -v $SPECIES
```

And finally, annotate the variants:
```
# 1. Further filter out intergenic variants, deal with alternative transcripts:
java -Xmx4g -jar ~/.conda/envs/SnpEff/SnpEff/snpEff.jar \
-c ~/.conda/envs/SnpEff/SnpEff/snpEff.config \
-v -no-downstream -no-intron -no-upstream -no-utr -no-intergenic \
-stats ${shortName}.filterstats.html \
-canon $SPECIES \
${base}.gatk.g.filt.vcf.gz > ${base}.gatk.g.filt.canon.ann.vcf && \
# 2. Extract the variant information required
java -Xmx4g -jar ~/.conda/envs/SnpEff/SnpEff/SnpSift.jar \
extractFields ${base}.gatk.g.filt.canon.ann.vcf \
CHROM ANN[0].GENEID POS REF ALT isHom"(GEN[0])" ANN[*].EFFECT \
> ${shortName}.just_genes_snps.txt
```

This will produce a file which ends in *.filterstats.genes.txt. This includes all types of variants, but the ones I want are variants_effect_missense_variant (nonsynonymous), variants_effect_stop_retained_variant and variants_effect_synonymous_variant (synonymous). There's also stuff like variants_effect_stop_gained or variants_effect_start_lost, which are clearly non-synonymous but we don't treat as non-synonymous variants because they're pseudogenising and may respond to selection differently compared to other non-synonymous changes. Similarly, we're not really interested in things like variants_effect_frameshift_variant or variants_effect_disruptive_inframe_deletion if we're calculating _dN/dS_.

I cut out the relevant columns with something like (column numbers may vary depending on how many different types of variants were annotated):
```
cat Mdes.filterstats.genes.txt | cut -f2,3,19,27,28 > Mdes.filterstats.genes.relevant.txt
```

### 4.4. Calculating _dNdS_

Now we have the ingredients required: counts of (i) synonymous and nonsynonymous sites and (ii) synonymous and nonsynonymous mutations.

I also produced the following files using the orthofinder outputs and the gifs:
- list of single-copy orthologs between species pairs
- which scaffolds correspond A and X for each focal species
- which genes are on which scaffold for each focal species

In RStudio, I combine all this to calcualte _dNdS_ on the X versus autosomes. I plot these with ggplot2 and run Mann-Whitney U tests to calculate significance. See code below.

```
### sex linkage/SCOG info

# load in SCOs for the focal species
SCOs <- read.table('single_copy_orthologs/bcop_all_SCO_isoforms.tsv', header=F)
colnames(SCOs) <- c('gene')
# sex linkage assignments for genes
linkage <- read.table('sex_linkage/assignments/bcop_sexlinkage.txt')
colnames(linkage) <- c('ctg', 'linkage')
# contig-gene info
ctg_chrom <- read.table('genes/Bcop.genes.tsv')
colnames(ctg_chrom) <- c('ctg', 'gene')
# merge ctg-gene with ctg-linkage
linkage <- merge(linkage, ctg_chrom, by=c('ctg'))
# now subset the genes by SCOs
linkage <- merge(linkage, SCOs, by=c('gene'))
nrow(linkage)
head(linkage)

### syn + nonsyn sites per gene

# load in degen counts
site_counts <- read.table('site_counts/bcop_degen_site_counts.tsv', header=T, stringsAsFactors=F)
site_counts <- site_counts[,c(1,6,7)]
colnames(site_counts) <- c('gene', 'N_nonsyn_sites', 'N_syn_sites')
nrow(site_counts)
head(site_counts)

# read in counts of non + syn variants per gene
var_counts <- read.table('variant_counts/Mdes.filterstats.genes.relevant.txt', header=T, stringsAsFactors=F)
var_counts <- var_counts[,c(2,3:5)]

var_counts <- var_counts %>% 
  group_by(TranscriptId) %>% 
  summarise(sum(variants_effect_missense_variant), sum(variants_effect_stop_retained_variant), sum(variants_effect_synonymous_variant))
var_counts <- as.data.frame(var_counts)
colnames(var_counts) <- c('TranscriptId', 'variants_effect_missense_variant', 'variants_effect_stop_retained_variant', 'variants_effect_synonymous_variant')

var_counts$N_syn_var <- var_counts$variants_effect_stop_retained_variant+var_counts$variants_effect_synonymous_variant
var_counts$N_nonsyn_var <- var_counts$variants_effect_missense_variant
var_counts <- var_counts[,c(1,5,6)]
colnames(var_counts) <- c('gene', 'N_syn_var', 'N_nonsyn_var')
# merge
degen_counts_linkage <- merge(degen_counts, linkage, by=c('gene'))
degen_counts_linkage <- degen_counts_linkage[,c(1,2,3,5)]
colnames(degen_counts_linkage) <- c('gene', 'N_nonsyn_sites', 'N_syn_sites', 'linkage')
syn_sites_and_var <- merge(var_counts, degen_counts_linkage, by=c('gene'))
# calc dN/dS per gene
# dN/dS is the ratio of the number of nonsynonymous substitutions per non-synonymous site to the number of synonymous substitutions per synonymous site
# i.e. (N_nonsyn_var/N_nonsyn_sites)/(N_syn_var/N_syn_sites)
syn_sites_and_var$dN <- syn_sites_and_var$N_nonsyn_var/syn_sites_and_var$N_nonsyn_sites
syn_sites_and_var$dS <- syn_sites_and_var$N_syn_var/syn_sites_and_var$N_syn_sites
syn_sites_and_var$dNdS <- syn_sites_and_var$dN/syn_sites_and_var$dS

syn_sites_and_var <- syn_sites_and_var[complete.cases(syn_sites_and_var$dNdS),]
syn_sites_and_var <- syn_sites_and_var[which(syn_sites_and_var$linkage=="A" | syn_sites_and_var$linkage=="X"),]
dndsplot <- ggplot(syn_sites_and_var, aes(x=linkage, y=dNdS, fill=linkage)) + 
  geom_boxplot(position=position_dodge(width=0.7), width=0.5, notch=T, outlier.shape = NA) +
  xlab("Chromosome") + ylab(expression(paste(italic("dN/dS")))) +
  scale_y_continuous(breaks = seq(0, 0.4, 0.2)) + 
  coord_cartesian(ylim = c(0, 0.4)) +
  scale_fill_manual("", breaks = c("A", "X"),
                    labels = c("Autosomes", "X"),
                    values = c("gray", "skyblue3")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.line = element_line(linewidth = 1),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color = NA),
        legend.background = element_rect(fill = 'transparent'),
        legend.box.background = element_rect(fill = 'transparent'))
print(dndsplot)
ggsave(file="Bcop_dNdS_ylim.svg", plot=dndsplot, width=2, height=2.8, units = "in")

# medians, standard deviation, etc.

# A median dN
median(syn_sites_and_var[which(syn_sites_and_var$linkage=="A"),][,c(8)])
# X median dN
median(syn_sites_and_var[which(syn_sites_and_var$linkage=="X"),][,c(8)])
# A median dS
median(syn_sites_and_var[which(syn_sites_and_var$linkage=="A"),][,c(9)])
# X median dS
median(syn_sites_and_var[which(syn_sites_and_var$linkage=="X"),][,c(9)])

# remove infinite values before calculating dN/dS medians (caused by genes where dS = 0 and dN > 0. These cases are rare)
syn_sites_and_var <- syn_sites_and_var[which(syn_sites_and_var$dNdS != Inf),]

# A median dNdS
median(syn_sites_and_var[which(syn_sites_and_var$linkage=="A"),][,c(10)]) # 
# X median dNdS
median(syn_sites_and_var[which(syn_sites_and_var$linkage=="X"),][,c(10)]) # 
# N A genes
nrow(syn_sites_and_var[which(syn_sites_and_var$linkage=="A"),]) 
# N X genes
nrow(syn_sites_and_var[which(syn_sites_and_var$linkage=="X"),]) 

# stat test
# dNdS Wilcoxin W
test <- syn_sites_and_var[,c(5,10)]
test_long <- test %>% group_by(linkage) %>% mutate(rn=row_number()) %>% spread(linkage, dNdS) %>% select(-rn)
test_long <- as.data.frame(test_long)
wilcox.test(test_long$A , test_long$X)


# A dNdS SD
sd(sdtest[which(sdtest$linkage=="A"),][,c(10)]) # 
# X dNdS SD
sd(sdtest[which(sdtest$linkage=="X"),][,c(10)]) # 

```












