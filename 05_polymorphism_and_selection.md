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

In RStudio, _pN/pS_ was then calculated, plotted, and tested in the same way that _dN/dS_ was calculated. Alpha was calculated as:

```
df$alpha <- 1 - ((df$dS*df$pN) / (df$pS+df$dS))/((df$pS*df$dN)/(df$pS+df$dS))
```

To calculate _DoS_ (Direction of Selection), I summed _dN_, _dS_, _pN_ and _pS_ for genes on the X and autosomes before calculations to increase power to estimate selection. I then used a permutation test to calculate significance and confident intervals (see below), as in Mongue et al. 2019; 2022; Mongue & Baird 2024.

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
# This is DoS, as per Stoletzki + Eyre-Walker
DoSCalc<-function(dN,dS,pN,pS)
{
  DoS <- sum(dN) / sum(dN+dS) - sum(pN) / sum(pN+pS)
  return(DoS)
}

# subset for each species + category of interest, in my case, linkage in the genome
poptest <- pop[which(pop$species=="ling"),]
popA<-poptest[which(poptest$chr=="A"),]
popX<-poptest[which(poptest$chr=="X"),]

# generate point estimates
Aa<-DoSCalc(popA$dN,popA$dS,popA$pN,popA$pS)
Xa<-DoSCalc(popX$dN,popX$dS,popX$pN,popX$pS)

#for the record
# Bcop:
#> Aa
#[1] 0.003780074
#> Xa
#[1] -0.01598161
# Ling:
#> Aa
#[1] -0.03455932
#> Xa
#[1] -0.0356238

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
  bsstata[i]<-abs((DoSCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(DoSCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
}
# sort all the permutation difference-of-alpha-values
bsdist<-sort(bsstata)
# p-value = the proportion of time we see a result as extreme or more so than the observed result based on chance alone.
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val 
# Bcop A v X: -Inf
# Ling A v X: 0.9042

#•#•#•# 95% confidence intervals

# we generate confidence intervals by taking the set of genes used to calculate the actual alpha point estimate, and then generating new alphas by sampling that WITH replacement. The idea is that our observations (genes) come from some distribution pool of values, and while we don't know the full range of that space, we do know that all our observations fall within it. So when we sample with replacement we might pick gene 1 twice and not sample gene 3 at all, etc. In effect this gives us a statement about the variability within the data (in stats term this is non-parametric bootstrapping)

# set a size
size<-10000
# initialise a vector for values
bsstata<-rep(0,size)
for(i in 1:size)
{
  # get CIs... replace=T is the key feature here
  bsap<-popA[sample(nrow(popA),size=nrow(popA),replace=T),]
  bsstata[i]<-DoSCalc(bsap$dN,bsap$dS,bsap$pN,bsap$pS)
}
# sort
bsdist<-sort(bsstata)
# bracket my confidence intervals around the point estimate just for sanity check
# take note here: if we want a 95% CI, we actually want to look at the values in the 2.5th percentile and the 97.5th percentile so that 5% of the values total fall outside it (if we did from 5% to 95% that would be a 90%CI)
bsspecs<-c(bsdist[250],Aa,bsdist[9750])
bsspecs
# Bcop A: 0.001192997 0.003780074 0.006360778
# Bcop X: -0.023270186 -0.015981606 -0.009017158
# Ling A: -0.04394752 -0.03455932 -0.02529753
# Ling X: -0.0511882 -0.0356238 -0.0205530

bsspecs_combined <- data.frame(species = c("bcop", "bcop", "ling", "ling"),
                               chromosome = c("A", "X", "A", "X"),
                               estimate = c(0.003780074, -0.015981606, -0.03455932, -0.0356238),
                               low_percentile = c(0.001192997, -0.023270186, -0.0511882, -0.0511882),
                               high_percentile = c(0.006360778, -0.009017158, -0.0205530, -0.0205530))
# plot:
ggplot(bsspecs_combined, aes(x=species, y=estimate, colour=chromosome)) +         
  geom_point(aes(colour = chromosome), position = position_dodge(0.75), size=4) +
  geom_errorbar(aes(ymin = low_percentile, ymax = high_percentile, colour = chromosome), position = position_dodge(0.75), width=0.5, linewidth = 1)+
  xlab("") + ylab(expression(paste(italic("DoS"))))+
  ylim(-0.06,0.06)+
  geom_hline(yintercept=c(0), linetype="dotted", size=1)+
  scale_x_discrete(labels=c("bcop" = expression(paste(italic("B. coprophila"))),
                            "ling" = expression(paste(italic("L. ingenua")))))+
  scale_colour_manual("",breaks=c('A', 'X'),
                    labels=c('A', 'X'),
                    values=c('gray', 'skyblue3'))+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.line=element_line(linewidth=1),
        axis.text.y=element_text(size=14, colour="black"),
        axis.text.x=element_text(size=12, colour="black"),
        axis.title=element_text(size=16))#,
```

