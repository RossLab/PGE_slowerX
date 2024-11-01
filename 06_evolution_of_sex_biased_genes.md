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
library(readr)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(PMCMRplus)
library(rstatix)
#•# ----------------------- FUNCTIONS

### quantile normalization function
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

### Function to read and merge count files from a specified directory
merge_count_files <- function(directory) {
  # Get the list of files in the directory
  file_list <- list.files(path = directory, pattern = "*.counts.tsv", full.names = TRUE)
  
  # Initialize an empty data frame
  merged_data <- NULL
  
  # Loop through each file and merge
  for (file in file_list) {
    # Read the current file
    current_data <- read_tsv(file, col_names = TRUE)
    
    if (is.null(merged_data)) {
      # If merged_data is NULL, initialize it with the first file's data
      merged_data <- current_data
    } else {
      # Merge the current data with the existing merged data
      merged_data <- full_join(merged_data, current_data, by = "target_id")
    }
  }
  
  return(merged_data)
}
#•# ----------------------- FILE PREP

### load in scaffold-chrom assignments and gene positions; assign genes
scaf_assignments <- read.table('../sex_linkage/assignments/bcop_sexlinkage.txt', header=F)
colnames(scaf_assignments) <- c('scaffold', 'chrom')
gene_positions <- read.table('../genes/Bcop.genes.tsv', header=F)
colnames(gene_positions) <- c('scaffold', 'gene')
gene_ids <- merge(gene_positions, scaf_assignments, by=c('scaffold'))
gene_ids <- gene_ids[gene_ids$chrom=="X" | gene_ids$chrom=="A",]
gene_ids <- gene_ids[,c(2,3)]

### make counts table
# Directory containing the count files
merged_data <- merge_count_files("Bcop/counts/")
merged_data <- as.data.frame(merged_data)

# rename columns:
colnames(merged_data) <- c('gene', 'm_g_1', 'f_g_1', 'm_g_2', 'f_g_2', 'f_g_3', 'f_r_1', 'f_r_2', 'f_r_3', 'f_s_1', 'f_s_2', 'f_s_3', 'm_g_3', 'm_r_1', 'm_r_2', 'm_r_3', 'm_s_1', 'm_s_2', 'm_s_3')
head(merged_data)

body <- merged_data[,c(1,17:19,10:12)]
gonad <- merged_data[,c(1,2,4,13,3,5,6)]
repro <- merged_data[,c(1,14:16,7:9)]

# write out counts tables
write.table(body, file="Bcop/counts/bcop_body_counts.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = F, col.names = T)
write.table(gonad, file="Bcop/counts/bcop_gonad_counts.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = F, col.names = T)
write.table(repro, file="Bcop/counts/bcop_repro_counts.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = F, col.names = T)

### PCA on all samples
library(DESeq2)

####
counts <- read.table("../../thesis_chapters/chapter_5/b_cop_DGE/counts_bcop_all.txt", header=T, stringsAsFactors=F)
colnames(counts) <- c('gene','MG1','MG2','MG3',
                          'MR1','MR2','MR3',
                          'MB1', 'MB2', 'MB3',
                          'FG1','FG2','FG3',
                          'FR1','FR2','FR3',
                          'FB1','FB2','FB3')
####
rownames(counts) <- counts$gene
counts <- counts[,c(2:19)]
counts <- counts[complete.cases(counts),]
count_matrix<-as.matrix(counts)
count_matrix<-round(count_matrix)
coldata <- read.table('DGE/coldata.txt', header=T)
dds = DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~ sex + tissue)
dds = dds[ rowMeans(counts(dds)) > 4, ] 
nrow(dds) # 15848/20498
rld = rlog(dds, blind=FALSE)
# PCA
counts = plotPCA(rld, intgroup = c("sex"), returnData=TRUE)
percentVar = round(100 * attr(counts, "percentVar"))
p <- ggplot(counts, aes(PC1, PC2, color=sex)) + geom_point(size=5) +
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
print(p)

```

Then measuring sex bias. I used the specificity metric, which basically describes what proportion of a gene's expression comes from males or females - 1 is 100% female, 0 is 100% male, 0.5 is equal expression between the sexes. I did this separately for the somatic gonads/reproductive tissue and the carcass, and then everything together. E.g.:

# quantile normalise counts
counts <- read.delim('Bcop/counts/bcop_body_counts.txt')
counts <- counts[complete.cases(counts),]
# make gene row name
row.names(counts) <- counts$gene
counts <- counts[,c(2:ncol(counts))]
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
nrow(counts_norm_means) # 20498
counts_norm_means <- counts_norm_means[rowSums(counts_norm_means[,c(2:3)] > 4) >= 1, ]
nrow(counts_norm_means) # 11597 body / 10506 repro
# calculate SPMs
count_SPMs <- counts_norm_means
count_SPMs$SPM <- (count_SPMs$female_mean^2)/((count_SPMs$male_mean^2)+(count_SPMs$female_mean^2))
# add linkage
nrow(count_SPMs) # 11597/10506
count_SPMs <- merge(count_SPMs, gene_ids, by=c('gene'))
nrow(count_SPMs) # 11597/10506

# assign biases
count_SPMs$sexbias <- "unbiased"
count_SPMs <- within(count_SPMs, sexbias[SPM > 0.7] <- 'female-biased')
count_SPMs <- within(count_SPMs, sexbias[SPM < 0.3] <- 'male-biased')
count_SPMs <- within(count_SPMs, sexbias[SPM > 0.9] <- 'strongly-female-biased')
count_SPMs <- within(count_SPMs, sexbias[SPM < 0.1] <- 'strongly-male-biased')

test <- count_SPMs[,c(1,5,6)]
table <- table(test$chrom, test$sexbias)
table
chisq.test(test$chrom, test$sexbias, correct=FALSE)
test<-chisq.test(test$chrom, test$sexbias, correct=FALSE)
test$residuals

# plot proportions
# Get the total genes per chromosome
gene_and_chrom_id <- counts_sexbias[,c(1,5)]
head(gene_and_chrom_id)
total_genes_per_chr <- as.data.frame(table(gene_and_chrom_id$chrom))
colnames(total_genes_per_chr) <- c("chrom", "total_genes")
# Make proportion of those genes which fall into a certain category
head(counts_sexbias)
bias_and_chrom_id <- counts_sexbias[,c(5,6)]
total_genes_diff_exp_by_chr <- as.data.frame(table(bias_and_chrom_id))
diff_exp_all_data <- merge(total_genes_diff_exp_by_chr, total_genes_per_chr, by ="chrom")
diff_exp_all_data$proportion <- diff_exp_all_data$Freq / diff_exp_all_data$total_genes
head(diff_exp_all_data)
# histogram showing proportion of biased genes per chromosome
class_order <- c("strongly-male-biased", "male-biased", "unbiased", "female-biased", "strongly-female-biased")
ggplot(diff_exp_all_data, aes(x=chrom, fill=factor(sexbias,level=class_order), y=proportion))+
  geom_bar(position = "stack", stat = "identity", colour="black")+
  ylab("Proportion of Genes")+
  xlab("Chromosome")+
  scale_fill_manual("",breaks=c("strongly-male-biased", "male-biased", "unbiased", "female-biased", "strongly-female-biased"),
                    labels=c("strongly-male-biased", "male-biased", "unbiased", "female-biased", "strongly-female-biased"),
                    values=c("darkblue", 'blue', 'grey', 'orange', "darkorange"))+
  scale_y_continuous(n.breaks = 3)+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.line=element_line(linewidth=1),
        axis.text=element_text(size=14, colour="black"),
        axis.title=element_blank())
```

I classified genes into 'predominantly expressed in gonads' vs 'predominantly expressed in body' to look at whether sexually dimorphic tissues were evolving differently.


## 5.2. Evolution of sex-biased genes

Here I essentially just used the assignments generated above, where I categorised genes as male-biased, unbiased etc. and combined this with the _dN/dS_ and _pN/pS_ results to look at evolution of differentially expressed genes.









