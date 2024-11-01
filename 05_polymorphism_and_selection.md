# 5. Polymorphism and selection

I used exactly the same pipeline for calling and annotating variants for the polymorphism data as for the divergence data - see 04_divergence.md. The only difference was that I combined VCFs from each resequening replicate prior to genotyping and filtering variants, i.e.:

```
# combine
/ceph/users/rbaird/software/gatk/gatk CombineGVCFs \
   -R $GENOME \
   --variant Bcop_reseq4.gatk.vcf.gz \
   --variant Bcop_reseq7.gatk.vcf.gz \
   --variant Bcop_reseq11.gatk.vcf.gz \
   -O Bcop.XX.reseq.cohort.vcf.gz

 # genotype
/ceph/users/rbaird/software/gatk/gatk --java-options "-Xmx4g" GenotypeGVCFs \
	-R $GENOME \
	-V Bcop.XX.reseq.cohort.vcf.gz \
	-O Bcop.XX.reseq.cohort.g.vcf.gz

# filter
/ceph/users/rbaird/software/gatk/gatk SelectVariants \
	-R $GENOME \
	-V Bcop.XX.reseq.cohort.g.vcf.gz \
	-O Bcop.XX.reseq.cohort.g.filt.vcf.gz \
	--select-type-to-include SNP \
	-select "QD > 2.0" \
	-select "FS < 60.0" \
	-select "MQ > 40.0"
```

Then SNPs were annotated and filtered with SnpEff and SnpSift as before.

In RStudio, _pN/pS_ was then calculated, plotted, and tested in the same way that _dN/dS_ was calculated (see divergence markdown). Alpha was calculated as:

```
# remove rows where both Ds and Ps = 0
df <- df[which((df$dS + df$pS)!=0),]
# calculate alpha
df$alpha <- 1 - ((df$dS*df$pN) / (df$pS+df$dS))/((df$pS*df$dN)/(df$pS+df$dS))

# plot:
p <- ggplot(data=pop, aes(x=species, y=alpha, fill=chr))+
  geom_boxplot(position=position_dodge(width=0.7), width=0.5, notch=T, outlier.size=0.5, outlier.alpha=0.5)+ 
  scale_y_continuous(breaks = seq(-5,1,1))  +
  coord_cartesian(ylim = c(-5, 1)) +
  theme_linedraw()+
  xlab("") + ylab("α (per-gene)")+
  scale_x_discrete(labels=c("bcop" = expression(paste(italic("B. coprophila"))),
                            "ling" = expression(paste(italic("L. ingenua"))),
                            "mdes" = expression(paste(italic("M. destructor")))))+
  scale_fill_manual("",breaks=c('A', 'X'),
                    labels=c('Autosomes', 'X'),
                    values=c('gray', 'skyblue3'))+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.line=element_line(linewidth=1),
        axis.text.y=element_text(size=12, colour="black"),
        axis.text.x=element_text(size=12, colour="black"),
        axis.title=element_text(size=16))#,
print(p)

# calculate medians:
pop <- pop[complete.cases(pop$alpha),]
pop <- pop[which(pop$alpha != -Inf),]
median(pop[which(pop$species=="bcop" & pop$chr=="A"),][,c('alpha')]) # 0.1632373
median(pop[which(pop$species=="bcop" & pop$chr=="X"),][,c('alpha')]) # 0.1428571
median(pop[which(pop$species=="ling" & pop$chr=="A"),][,c('alpha')]) # 0.0939976
median(pop[which(pop$species=="ling" & pop$chr=="X"),][,c('alpha')]) # 0.2496499
median(pop[which(pop$species=="mdes" & pop$chr=="A"),][,c('alpha')]) # 0.2381294
median(pop[which(pop$species=="mdes" & pop$chr=="X"),][,c('alpha')]) # 0.3333333

# stat test
test <- pop[which(pop$species=="mdes"),]
test <- test[,c(3,8)]
# spread to make partitions into columns
test_long <- test %>% group_by(chr) %>% mutate(rn=row_number()) %>% spread(chr, alpha) %>% select(-rn)
test_long <- as.data.frame(test_long)
colnames(test_long) <- c('A', 'X')

wilcox.test(test_long$A , test_long$X)
```

To calculate aggregated alpha, I summed _dN_, _dS_, _pN_ and _pS_ for genes on the X and autosomes before calculations to increase power to estimate selection. I then used a permutation test to calculate significance and confident intervals (see below), as in Mongue et al. 2019; 2022; Mongue & Baird 2024.

```
# I first created a dataframe with _dN_, _dS_, _pN_ and _pS_ for the two species being tested.
pop <- merge(bpop, lpop, by=c('species', 'chr', 'dN', 'dS', 'pN', 'pS'), all=T)

# remove rows where both Ds and Ps = 0
pop<-pop[which((pop$dS + pop$pS)!=0),]

########## Point estimates & permutations
#define a function for calculating alpha/NI/DoS
#note this calculates the Tarone/Greenland neutrality index, so alpha is 1 - NItgCalc in this case
NItgCalc<-function(dN,dS,pN,pS)
{
  unbiasedNI<-sum((dS*pN)/(pS+dS))/sum((pS*dN)/(pS+dS))
  return(unbiasedNI)
}

# subset for each species + category of interest, in my case, linkage in the genome
poptest <- pop[which(pop$species=="ling"),]
popA<-poptest[which(poptest$chr=="A"),]
popX<-poptest[which(poptest$chr=="X"),]

# generate point estimates
Aa<-1-NItgCalc(popA$dN,popA$dS,popA$pN,popA$pS)
Xa<-1-NItgCalc(popX$dN,popX$dS,popX$pN,popX$pS)

#for the record
# Bcop:
#> Aa
#[1] 0.05861344
#> Xa
#[1] -0.06999781
# Ling:
#> Aa
#[1] -0.1573217
#> Xa
#[1] -0.1252326
# Mdes:
#> Aa
#[1] 0.1426937
#> Xa
#[1] 0.1478854


#•#•#•# permute

# test statistic = the absolute value of the difference of the two point estimates
testor<-abs((Aa-Xa))
# set a number of permutations
size<-10000
# pre-initialise a vector to store all our test statistics (i.e. the difference between random draws)
bsstata<-rep(0,size)
# make the pool of genes to draw from by row-binding the two categories
mwhole<-rbind(popA,popX)
# pull out two groups equal to the two sizes from the real statistic:
pa<-nrow(popA)
pb<-nrow(mwhole)
# run the permutation
for(i in 1:size)
{
  # shuffle our combined gene set, without reusing genes
  bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  # then we do the difference of alpha calc: what's happening here is we're calculating the absolute difference in alpha for two groups of size pa and pb, as defined above. Note that you don't actually have to split the "gene pool" into two objects, you can use row indexing to just do the calcs on the relevant parts
  bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
}
# sort all the permutation difference-of-alpha-values
bsdist<-sort(bsstata)
# p-value = the proportion of time we see a result as extreme or more so than the observed result based on chance alone.
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val 
# Bcop A v X: -Inf
# Ling A v X: 0.5262
# Mdes A v X: 0.7782

#•#•#•# 95% confidence intervals

# we generate confidence intervals by taking the set of genes used to calculate the actual alpha point estimate, and then generating new alphas by sampling that WITH replacement. The idea is that our observations (genes) come from some distribution pool of values, and while we don't know the full range of that space, we do know that all our observations fall within it. So when we sample with replacement we might pick gene 1 twice and not sample gene 3 at all, etc. In effect this gives us a statement about the variability within the data (in stats term this is non-parametric bootstrapping)

# set a size
size<-10000
# initialise a vector for values
bsstata<-rep(0,size)
for(i in 1:size)
{
  # get CIs... replace=T is the key feature here
  bsap<-popX[sample(nrow(popX),size=nrow(popX),replace=T),]
  bsstata[i]<-1-NItgCalc(bsap$dN,bsap$dS,bsap$pN,bsap$pS)
}
# sort
bsdist<-sort(bsstata)
# bracket my confidence intervals around the point estimate just for sanity check
# take note here: if we want a 95% CI, we actually want to look at the values in the 2.5th percentile and the 97.5th percentile so that 5% of the values total fall outside it (if we did from 5% to 95% that would be a 90%CI)
bsspecs<-c(bsdist[250],Xa,bsdist[9750])
bsspecs
# Bcop A: 0.03818344 0.05861344 0.07849805
# Bcop X: -0.12608776 -0.06999781 -0.01771537
# Ling A: -0.2009956 -0.1573217 -0.1144948
# Ling X: -0.20873469 -0.12523260 -0.04685418
# Mdes A: 0.1180767 0.1426937 0.1663068
# Mdes X: 0.1200380 0.1478854 0.1744304

bsspecs_combined <- data.frame(species = c("bcop", "bcop", "ling", "ling", "mdes", "mdes"),
                               chromosome = c("A", "X", "A", "X", "A", "X"),
                               estimate = c(0.05861344, -0.06999781, -0.1573217, -0.12523260, 0.1426937, 0.1478854),
                               low_percentile = c(0.03818344, -0.12608776, -0.2009956, -0.20873469, 0.1180767, 0.1200380),
                               high_percentile = c(0.07849805, -0.01771537, -0.1144948, -0.04685418, 0.1663068, 0.1744304))
# plot:
p <- ggplot(bsspecs_combined, aes(x=species, y=estimate, colour=chromosome)) +         
  geom_point(aes(colour = chromosome), position = position_dodge(0.75), size=4) +
  geom_errorbar(aes(ymin = low_percentile, ymax = high_percentile, colour = chromosome), position = position_dodge(0.75), width=0.5, linewidth = 1)+
  xlab("") + ylab("α (aggregate)") +
  scale_y_continuous(limits = c(-0.25, 0.25), breaks = seq(-0.25,0.25,0.1))  +
  geom_hline(yintercept=c(0), linetype="dotted", size=1)+
  scale_x_discrete(labels=c("bcop" = expression(paste(italic("B. coprophila"))),
                            "ling" = expression(paste(italic("L. ingenua"))),
                            "mdes" = expression(paste(italic("M. destructor")))))+
  scale_colour_manual("",breaks=c('A', 'X'),
                      labels=c('Autosomes', 'X'),
                      values=c('gray', 'skyblue3'))+
  theme_classic()+
  #theme(legend.position = "none")+
  theme(axis.line=element_line(linewidth=1),
        axis.text.y=element_text(size=14, colour="black"),
        axis.text.x=element_text(size=12, colour="black"),
        axis.title=element_text(size=16),
        legend.text=element_text(size=12))#,
print(p)
```

