# 5. Evolution of sex-biased genes

## 5.1. Differential gene expression

First I want to look at sex-biased gene content on the sex chromosomes, as this may have an effect on how they evolve. For example, if there are lots of male-biased genes then we may expect the contribution of haploid selection to be stronger.

I have the following RNAseq datasets:

| Species | Sex | Tissue | Replicates
| - | - | - | - |
| Bradysia coprophila | Male | Testes | 3 |
| Bradysia coprophila | Male | Ejaculatory Sac | 3 |
| Bradysia coprophila | Male | Carcass | 3 |
| Bradysia coprophila | Androgenic female | Ovaries | 3 |
| Bradysia coprophila | Androgenic female | Reproductive glands | 3 |
| Bradysia coprophila | Androgenic female | Carcass | 3 |
| Bradysia coprophila | Gynogenic female | Ovaries | 3 |
| Bradysia coprophila | Gynogenic female | Reproductive glands | 3 |
| Bradysia coprophila | Gynogenic female | Carcass | 3 |
| Lycoriella ingenua | Male | Testes | 3 |
| Lycoriella ingenua | Male | Ejaculatory Sac | 3 |
| Lycoriella ingenua | Male | Carcass | 3 |
| Lycoriella ingenua | Female | Ovaries | 3 |
| Lycoriella ingenua | Female | Reproductive glands | 3 |
| Lycoriella ingenua | Female | Carcass | 3 |

The reason I separated reproductive tissues into germline and non-germline is because the germline is XX, so presumably won't drive faster-X evolution.

Gynogenic (female-producer) Bradysia coprophila females are X'X, i.e. they have a huge female-limited heterozygous inversion. This will affect measured of sex-biased gene expression, so I will only compare androgenic (male-producer, XX) females with males. Lycoriella ingenua doesn't have this problem - all females are XX.

To quantify gene expression I used Kallisto. I pulled out the longest isoforms from the predicted transcripts from the Bradysia coprophila and Lycoriella annotations:

```
cat augustus_hints.fasta | awk '/^>/ {if(N>0) printf("\n"); printf("%s\t",$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' | tr "." "\t" | awk -F '	'	'{printf("%s\t%d\n",$0,length($3));}' | sort -t '	' -k1,1 -k4,4nr | sort -t '	' -k1,1 -u -s | sed 's/	/./' | cut -f 1,2 | tr "\t" "\n"	| fold -w 60 > augustus.longest.fasta
```

And then ran kallisto to get counts per transcript:

```
kallisto index -i longest_trx.idx augustus.longest.fasta

for file in $(ls *_1.trimmed.fq.gz)
do
	base=$(basename ${file} "_1.trimmed.fq.gz")
	kallisto quant \
	--index=longest_trx.idx \
	-b 20 \
	--output-dir=${base}_out \
	--threads=16 \
	--plaintext ${base}_1.trimmed.fq.gz ${base}_2.trimmed.fq.gz \
	&& mv ${base}_out/abundance.tsv ${base}.tsv
done
```

In RStudio I peformed a bit of QC on the datasets - I plotPCA from the DESeq2 package and also compared TPM distributions (figures are in supplementary materials). I also quantile-normalised the counts.

To make reading files into R easier, I ran a small bash script which will reformat the *.tsv files that were output by Kallisto. It will put the file name as a column and remove columns that won't be needed:

```
for file in $(ls *.tsv)
do
	base=$(basename ${file} ".tsv")
	cat ${base}.tsv \
	| cut -f1,4 \
	| sed "s/est_counts/${base}_counts/" \
	> ${base}.counts.mod.tsv
done
```

Now in RStudio, loading in the results, running the PCA, and defining the quantile normalisation function:

```
library(data.table)
library(tidyverse)
library(dplyr)

# Generate a counts table
all_files <- dir("/Users/robertbaird/documents/analyses/chapter_5/ling_DGE_2")
file_names <- grep(all_files,pattern = "*.counts.mod.tsv",value = TRUE)
Data_file <- map(file_names,read.delim, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
Merge_All_Samples <- Data_file %>% reduce(inner_join, by = "target_id")
head(Merge_All_Samples)
# L1 L10 L11 L12 L13 L14 L15 L16 L17 L18 L2 L3 L4 L5 L6 L7 L8 L9
Merge_All_Samples <- Merge_All_Samples[,c(1,2,12:19,3:11)]
colnames(Merge_All_Samples) <- c('gene','male_gonad_1','male_gonad_2','male_gonad_3',
                                 'male_repro_tract_1','male_repro_tract_2','male_repro_tract_3',
                                 'male_carcass_1', 'male_carcass_2', 'male_carcass_3',
                                 'female_gonad_1','female_gonad_2','female_gonad_3',
                                 'female_repro_tract_1','female_repro_tract_2','female_repro_tract_3',
                                 'female_carcass_1','female_carcass_2','female_carcass_3')

counts_ling_gonads <- Merge_All_Samples[,c(1:4,11:13)]
counts_ling_repro_tract <- Merge_All_Samples[,c(1,5:7,14:16)]
counts_ling_carcass <- Merge_All_Samples[,c(1,8:10,17:19)]
# write out counts tables
write.table(counts_ling_gonads, file="counts_ling_gonads.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = F, col.names = T)
write.table(counts_ling_repro_tract, file="counts_ling_repro_tract.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = F, col.names = T)
write.table(counts_ling_carcass, file="counts_ling_carcass.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = F, col.names = T)
write.table(Merge_All_Samples, file="counts_ling_all.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = F, col.names = T)

library(DESeq2)
library(apeglm)
library(reshape2)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(PoiClaClu)

### PCA on all samples

all_counts <- read.table("counts_ling_all.txt", header=T, stringsAsFactors=F)
colnames(all_counts) <- c('gene','MG1','MG2','MG3',
                          'MR1','MR2','MR3',
                          'MB1', 'MB2', 'MB3',
                          'FG1','FG2','FG3',
                          'FR1','FR2','FR3',
                          'FB1','FB2','FB3')
#
rownames(all_counts) <- all_counts$gene
all_counts <- all_counts[,c(2:19)]
all_counts <- all_counts[complete.cases(all_counts),]
count_matrix<-as.matrix(all_counts)
count_matrix<-round(count_matrix)
coldata <- read.table('coldata.txt', header=T)
dds = DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~ sex + tissue)
dds = dds[ rowMeans(counts(dds)) > 4, ] 
nrow(dds) # 22805/31796
rld = rlog(dds, blind=FALSE)
# PCA
counts = plotPCA(rld, intgroup = c("sex"), returnData=TRUE)
percentVar = round(100 * attr(counts, "percentVar"))
ggplot(counts, aes(PC1, PC2, color=sex)) + geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()+
  geom_text_repel(aes(label=c('MG1','MG2','MG3',
                              'MR1','MR2','MR3',
                              'MB1', 'MB2', 'MB3',
                              'FG1','FG2','FG3',
                              'FR1','FR2','FR3',
                              'FB1','FB2','FB3')), 
                  size=6,show.legend=FALSE, point.padding = 1.0, box.padding = 0.5, max.overlaps=100) +
  theme_bw()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18))

## quantile normalization function

quantile_normalisation <- function(df){
  # Determine the ranks of each column from lowest to highest
  df_rank <- apply(df,2,rank,ties.method="min")
  # Sort the original matrix from lowest to highest
  df_sorted <- data.frame(apply(df, 2, sort))
  # Calculate the means
  df_mean <- apply(df_sorted, 1, mean)
  # Substitute the means into the ranked matrix
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}
```

Then measuring sex bias. I used the specificity metric, which basically describes what proportion of a gene's expression comes from males or females - 1 is 100% female, 0 is 100% male, 0.5 is equal expression between the sexes. I did this separately for the somatic gonads/reproductive tissue and the carcass, and then everything together. E.g.:

```
### GONADS

counts_repro <- read.delim('ling_DGE/counts_ling_repro_tract.txt')
counts <- counts_repro[complete.cases(counts_repro),]
# make gene row name
row.names(counts) <- counts$gene
counts <- counts[,c(2:7)]
# quantile normalisation
counts_norm <- quantile_normalisation(counts)
counts_norm <- as.data.frame(counts_norm)
counts_norm$gene <-rownames(counts_norm)
# calculate means
counts_norm_means <- counts_norm
counts_norm_means$male_mean <- rowMeans(counts_norm_means[,c(1:3)])
counts_norm_means$female_mean <- rowMeans(counts_norm_means[,c(4:6)])
counts_norm_means <- counts_norm_means[,c(7:9)]
# normalised count of >4 in at least one tissue in at least one sex (i.e. filter out noise)
counts_norm_means <- counts_norm_means[rowSums(counts_norm_means[,c(2:3)] > 4) >= 1, ] # 
# calculate SPMs
count_SPMs <- counts_norm_means
count_SPMs$SPM <- (count_SPMs$female_mean^2)/((count_SPMs$male_mean^2)+(count_SPMs$female_mean^2))
# add linkage
count_SPMs <- merge(count_SPMs, gene_ids, by=c('gene')) # 
# assign biases
count_SPMs$sexbias <- "unbiased"
count_SPMs <- within(count_SPMs, sexbias[SPM > 0.7] <- 'female-biased')
count_SPMs <- within(count_SPMs, sexbias[SPM < 0.3] <- 'male-biased')
count_SPMs <- within(count_SPMs, sexbias[SPM > 0.9] <- 'strongly-female-biased')
count_SPMs <- within(count_SPMs, sexbias[SPM < 0.1] <- 'strongly-male-biased')

sexbias_sex_tissues <- count_SPMs
write.table(sexbias_sex_tissues, file="ling_DGE/ling_gonads_somatic_res_SPM.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = F, col.names = T)

# chi sq test
head(sexbias_sex_tissues)

test <- sexbias_sex_tissues[,c(1,5,6)]
table <- table(test$chrom, test$sexbias)
table
#female-biased male-biased strongly-female-biased strongly-male-biased unbiased
#A          1499        1683                   3290                 3357     2525
#X           377         368                    767                  832      579
chisq.test(test$chrom, test$sexbias, correct=FALSE)
# A and X do NOT differ in sex-biased composition  (X-squared = 5.0373, df = 4, p-value = 0.2835)
# check residuals
test<-chisq.test(test$chrom, test$sexbias, correct=FALSE)
test$residuals
#test$chrom female-biased male-biased strongly-female-biased strongly-male-biased   unbiased
#A    -0.4636418   0.5997395              0.1613101           -0.5241116  0.2973759
#X     0.9531734  -1.2329687             -0.3316279            1.0774897 -0.6113573

# plot proportions
# Get the total genes per chromosome
gene_and_chrom_id <- sexbias_sex_tissues[,c(1,5)]
head(gene_and_chrom_id)
total_genes_per_chr <- as.data.frame(table(gene_and_chrom_id$chrom))
colnames(total_genes_per_chr) <- c("chrom", "total_genes")
# Make proportion of those genes which fall into a certain category
head(sexbias_sex_tissues)
bias_and_chrom_id <- sexbias_sex_tissues[,c(5,6)]
total_genes_diff_exp_by_chr <- as.data.frame(table(bias_and_chrom_id))
diff_exp_all_data <- merge(total_genes_diff_exp_by_chr, total_genes_per_chr, by ="chrom")
diff_exp_all_data$proportion <- diff_exp_all_data$Freq / diff_exp_all_data$total_genes
head(diff_exp_all_data)
# histogram showing proportion of biased genes per chromosome
class_order <- c("male-limited", "male-biased", "unbiased", "female-biased", "female-limited")
ggplot(diff_exp_all_data, aes(x=chrom, fill=factor(sexbias,level=class_order), y=proportion))+
  geom_bar(position = "stack", stat = "identity", colour="black")+
  ylab("Proportion of Genes")+
  xlab("Chromosome")+
  scale_fill_manual("",breaks=c("male-limited", "male-biased", "unbiased", "female-biased", "female-limited"),
                    labels=c("male-limited", "male-biased", "unbiased", "female-biased", "female-limited"),
                    values=c("darkblue", 'blue', 'grey', 'orange', "darkorange"))+
  scale_y_continuous(n.breaks = 3)+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.line=element_line(linewidth=1),
        axis.text=element_text(size=14, colour="black"),
        axis.title=element_blank())
```


## 5.2. Evolution of sex-biased genes

Here I essentially just used the assignments generated above, where I categorised genes as male-biased, unbiased etc. and combined this with the _dN/dS_ and _pN/pS_ results to look at evolution of differentially expressed genes.









