# Sex linkage assignments

Where possible, I used male reads to assign sex-linked scaffolds. This wasn't necessary for Bradysia coprophila, for which there is already a chromosome-level assembly with an X chromosome scaffold (Urban et al. 2022). For all other species I needed to assign X-linked scaffolds. This was done via alignments of male reads to their respective genomes with BWA:

```
bwa mem -t 32 ${base}.*.fa ${base}_1.trimmed.fastq.gz ${base}_2.trimmed.fastq.gz | \
	samtools view -b - > ${base}.bwa.bam && \
	samtools sort -@ 8 -o ${base}.bwa.sorted.bam ${base}.bwa.bam
```

For one species, Aphidoletes aphidimyza, I used HiFi data generated for a single male at the Ross Lab, University of Edinburgh, to ID sex-linked scaffolds. This is the only sequencing library used for this study which is not being made publicly available as it is the focus of a future project in the Ross Lab. This alignment was produced with minimap2:

```
minimap2 -t 32 -ax map-pb $GENOME $READS | samtools view -bo $READS.bam && samtools sort -o $READS.sorted.bam $READS.bam
```

Then, I used samtools coverage to get read depth for each scaffold:

```
samtools coverage ${base}.bam | cut -f1,2,3,7 > ${base}.cov
```

The resulting coverage file was then analysed in RStudio to assign sex-linked scaffolds. For chromosome-level assemblies (Aphidoletes aphidimyza and Sitodiplosis mosellana), this was simple as it requried assigning two scaffolds as autosomal and two as X-linked in each case. For others, I plotted covergae histograms, identified the haploid (1n) and diploid (2n) peaks, and assigned all scaffolds that lay 25% to the left and right of each peak as X-linked and autosomal. For example, if the 1n peak was at 100x and the 2n peak was at 200x, then scaffolds with 70-125x coverage were assigned as X-linked and 150-250x as autosomal. Anything else was left unassigned. All histograms with the regions assigned highlighted, as well as the amount of sequence assigned in each species, are available in the supplementary materials. 


For the following species, only female reads were available: Bradysia odoriphaga, Lycoriella agraria, Mayetiola hordei, Contarinia nasturtii, Contarinia rumicis. Sex-linked scaffolds for these species (except Contarinia) were assigned using alignments to the closest outgroup (i.e. their species pair - see 04_divergence.md). This was done with minimap2, e.g..:

```
# align with minimap2
minimap2 -t 40 -c --secondary=no -ax asm5 species1.fa species2.fa | samtools view -bo sp1_v_sp2.bam && samtools sort -o sp1_v_sp2.sorted.bam sp1_v_sp2.bam && rm sp1_v_sp2.bam
# convert bam to bed
bedtools bamtobed -ed -i sp1_v_sp2.sorted.bam > sp1_v_sp2.sorted.bed && rm sp1_v_sp2.sorted.bam
```

Using the bed output, I assigned species2 scaffolds based on alignments to species1 scaffolds with >70% identity. If a species2 scaffold aligned to multiple species1 scaffolds then I took the longest alignment. This was all done with a some simple R code:

```
aln <- read.table('bodo_v_bcop.sorted.bed', header=F, stringsAsFactors=F)
colnames(aln) <- c('ctg_ref', 'aln_end', 'aln_start', 'ctg_query', 'N_mismatches', 'orientation')
aln$aln_len <- aln$aln_start-aln$aln_end
head(aln)

# sum the aln length for each contig
getting_len <- aln %>%
  group_by(ctg_ref, ctg_query) %>%
  summarise(total_aln_len = sum(aln_len))
head(getting_len)
nrow(getting_len)
# sum the n mismatches for each contig
getting_NM <- aln %>%
  group_by(ctg_ref, ctg_query) %>%
  summarise(total_mismatches = sum(N_mismatches))
head(getting_NM)
nrow(getting_NM)
# merge and calculate percentage identity for each contig vs the core genome
total_aln<- merge(getting_NM, getting_len, by=c('ctg_ref', 'ctg_query'))
head(total_aln)
total_aln$pident <- 1-total_aln$total_mismatches/total_aln$total_aln_len
# get the top hit for each scaffold
best_hits <- total_aln %>% 
  dplyr::group_by(ctg_query) %>% 
  dplyr::slice(which.max(total_aln_len))
head(best_hits)
nrow(best_hits)
best_hits <- best_hits[which(best_hits$pident > 0.70),]
nrow(best_hits)
# here I'm assigning scaffolds to chromosomes based on the longest alignment, and only with >70% identity

# outgroup  assignments
outgr_assn <- read.table('../assignments/bcop_sexlinkage.txt', header=F)
colnames(outgr_assn) <- c('ctg_ref', 'linkage')
best_hits <- merge(best_hits, outgr_assn, by=c('ctg_ref')) 
linkage <- best_hits[,c(2,6)]

# ctg lengths
lengths <- read.table('genome_lengths/Bodo.ncbi.fa.lengths', header=F, stringsAsFactors=F)
colnames(lengths) <- c('ctg_query', 'len')

linkage_lengths <- merge(linkage, lengths, by=c("ctg_query"))
sum(linkage_lengths[which(linkage_lengths$linkage=="A"),]$len) # 263514572bp
sum(linkage_lengths[which(linkage_lengths$linkage=="X"),]$len) # 78764338bp
```

For one species pair - Contarinia nasturtii and Contarinia rumicis, male reads were not available for either species. I assigned scaffolds from these species based on alignments to the closest outgroup - in this case Sitodiplosis mosellana.






