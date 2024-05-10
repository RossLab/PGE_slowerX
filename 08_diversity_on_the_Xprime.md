# 8. Diversity on the B coprophila X

I estimated dN/dS on the X' supergene by following the same approach as in the _dN/dS_ analysis for Bradysia coprophila, but mapping the outgroup (Bradysia odoriphaga) to the supergene scaffold. 

I also did a sliding window analysis of heterozygosity along the X with two of the population resequencing samples. To do this, I called variants with GATK-4 as in the _pN/pS_ analysis, but this time including the options `--output-mode EMIT_ALL_CONFIDENT_SITES` and `--include-non-variant-sites` so to call genotypes at every site (both variant and invariant):

```
# call variants
/ceph/users/rbaird/software/gatk/gatk --java-options "-Xmx4g" HaplotypeCaller \
-R $GENOME \
-I ${base}.pic.rmdup.rg.sorted.bam \
-O ${base}.gatk.vcf.gz \
--output-mode EMIT_ALL_CONFIDENT_SITES \
-ERC GVCF && \
# genotype variants
/ceph/users/rbaird/software/gatk/gatk --java-options "-Xmx4g" GenotypeGVCFs \
-R $GENOME \
-V ${base}.gatk.vcf.gz \
-O ${base}.gatk.g.vcf.gz \
--include-non-variant-sites
```

Then I used the popgenWindows.py script from the genomics_general toolkit to calcualte heterozygosity across 100kb windows:

```
python ~/software/genomics_general/VCF_processing/parseVCF.py -i ${base}.gatk.g.vcf.gz --skipIndels --minQual 30 --gtf flag=DP min=4 max=50 -o ${base}.snps.geno.gz
python ~/software/genomics_general/filterGenotypes.py --threads 8 -i ${base}.snps.geno.gz -o ${base}.snps.geno.filt.gz --minAlleles 1 --maxAlleles 2 -of alleles
python ~/software/genomics_general/popgenWindows.py -w 100000 --analysis indHet -g ${base}.snps.geno.gz -o ${base}indHetdens.100kbwins.csv.gz -f phased -T 8
```