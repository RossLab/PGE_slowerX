# Genome assembly and annotation

## Genome assembly

De novo short-read (Illumina) genome assemblies were generated for _Bradysia confinis_, _Bradysia desolata_, _Bradysia pectoralis_, _Lycoriella ingenua_ and _Lycoriella agraria_. This was done using SPAdes:

```
base=$(basename $file "_1.trimmed.fastq.gz")
spades.py -t 16 --isolate -1 ${base}_1.trimmed.fastq.gz -2 ${base}_2.trimmed.fastq.gz -o ${base}.spades
```

I then ran blobtools on the scaffolds: mapped the reads back to the assembly with bwa mem, ran blastn and diamond searches against the NCBI and Uniprot databases, and then generated the plots and hits tables with blobtools:

```
# 1. Mapping
bwa index ${base}.spades.scaffolds.fasta
bwa mem -t 32 ${base}.spades.scaffolds.fasta ${base}_1.trimmed.fastq.gz ${base}_2.trimmed.fastq.gz | \
samtools view -b - > ${base}.bwa.bam && \
samtools sort -@ 8 -o ${base}.bwa.sorted.bam ${base}.bwa.bam && \
samtools index ${base}.bwa.sorted.bam

# 2. Blastn
blastn -num_threads 30 -max_target_seqs 10 -max_hsps 1 -db /ceph/software/databases/ncbi/nt.56 -evalue 1e-25 -outfmt '6 qseqid staxids bitscore std' -query ${base}.spades.scaffolds.fasta -out ${base}.vs.ncbi.nt.out

# 3. Diamond 
diamond blastx --threads 30 --db /ceph/software/databases/uniprot_2018_04/full/reference_proteomes.dmnd --sensitive --max-target-seqs 1 --evalue 1e-25 --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --query ${base}.spades.scaffolds.fasta > ${base}.vs.uniprot.out

# 4. Blobtools
~/software/blobtools/blobtools create -i ${base}.spades.scaffolds.fasta -t ${base}.vs.ncbi_2019_09.nt.mts10.out -t ${base}.vs.uniprot.out -b ${base}.bwa.sorted.bam -o ${base}
~/software/blobtools/blobtools view -i ${base}.blobDB.json
~/software/blobtools/blobtools plot -i ${base}.blobDB.json
```

In RStudio, using the bestsum tables, I removed scaffolds with hits to non-metazoans as well as those with length <1kb and coverage <2x. Seqtk was then used to remove these scaffolds from the assemblies:

```
seqtk subseq ${base}.spades.scaffolds.fasta ${base}_contig_keep_list.txt > ${base}.spades.scaffolds.filt.fasta
```

Finally, I ran BUSCO v5 to assess assembly completeness:

```
busco -c 60 -i ${base}.fasta -l insecta_odb10 -o ${base}.busco -m genome
```

For _Lycoriella ingenua_, I sought to improve the assemly because I was running further downstream polymorphism and expression analyses with this species. I used the RNAseq data from this species to improve contiguity of coding regions. To do this, I first mapped the reads with TopHat:

```
bowtie2-build Lycoriella_ingenua_scaffolds.filt.fasta Lycoriella_ingenua_scaffolds.idx
# testes
tophat \
-p 16 \
-o testes.tophat_out \
Lycoriella_ingenua_scaffolds.idx \
L1_EKRN230017318-1A_HY5YHDSX5_L2_1.trimmed.fq.gz,L2_EKRN230017319-1A_HY5YHDSX5_L2_1.trimmed.fq.gz,L3_EKRN230017320-1A_HY5YHDSX5_L2_1.trimmed.fq.gz,L10_combined_1.trimmed.fq.gz,L11_combined_1.trimmed.fq.gz,L12_EKRN230017329-1A_HYCMYDSX5_L3_1.trimmed.fq.gz \
L1_EKRN230017318-1A_HY5YHDSX5_L2_2.trimmed.fq.gz,L2_EKRN230017319-1A_HY5YHDSX5_L2_2.trimmed.fq.gz,L3_EKRN230017320-1A_HY5YHDSX5_L2_2.trimmed.fq.gz,L10_combined_2.trimmed.fq.gz,L11_combined_2.trimmed.fq.gz,L12_EKRN230017329-1A_HYCMYDSX5_L3_2.trimmed.fq.gz
# ejaculatory sac
tophat \
-p 16 \
-o ejaculatory_sac.tophat_out \
Lycoriella_ingenua_scaffolds.idx \
L4_EKRN230017321-1A_HYCJ7DSX5_L4_1.trimmed.fq.gz,L5_EKRN230017322-1A_HY5YHDSX5_L2_1.trimmed.fq.gz,L6_EKRN230017323-1A_HY5YHDSX5_L2_1.trimmed.fq.gz,L13_EKRN230017330-1A_HYCMYDSX5_L3_1.trimmed.fq.gz,L14_EKRN230017331-1A_HYCMYDSX5_L3_1.trimmed.fq.gz,L15_combined_1.trimmed.fq.gz \
L4_EKRN230017321-1A_HYCJ7DSX5_L4_2.trimmed.fq.gz,L5_EKRN230017322-1A_HY5YHDSX5_L2_2.trimmed.fq.gz,L6_EKRN230017323-1A_HY5YHDSX5_L2_2.trimmed.fq.gz,L13_EKRN230017330-1A_HYCMYDSX5_L3_2.trimmed.fq.gz,L14_EKRN230017331-1A_HYCMYDSX5_L3_2.trimmed.fq.gz,L15_combined_2.trimmed.fq.gz \
# carcass
tophat \
-p 16 \
-o carcass.tophat_out \
Lycoriella_ingenua_scaffolds.idx \
L7_EKRN230017324-1A_HYCMYDSX5_L3_1.trimmed.fq.gz,L8_EKRN230017325-1A_HYCMYDSX5_L3_1.trimmed.fq.gz,L9_EKRN230017326-1A_HYCMYDSX5_L3_1.trimmed.fq.gz,L16_EKRN230017333-1A_HYCMYDSX5_L3_1.trimmed.fq.gz,L17_EKRN230017334-1A_HYCMYDSX5_L3_1.trimmed.fq.gz,L18_EKRN230017335-1A_HYCMYDSX5_L3_1.trimmed.fq.gz \
L7_EKRN230017324-1A_HYCMYDSX5_L3_2.trimmed.fq.gz,L8_EKRN230017325-1A_HYCMYDSX5_L3_2.trimmed.fq.gz,L9_EKRN230017326-1A_HYCMYDSX5_L3_2.trimmed.fq.gz,L16_EKRN230017333-1A_HYCMYDSX5_L3_2.trimmed.fq.gz,L17_EKRN230017334-1A_HYCMYDSX5_L3_2.trimmed.fq.gz,L18_EKRN230017335-1A_HYCMYDSX5_L3_2.trimmed.fq.gz \
```

Then used the alignments to improve scaffold with rascaf:
```
rascaf -b ${base}.bam -f Lycoriella_ingenua_scaffolds.filt.fasta -o ${base}
rascaf-join \
-r carcass_accepted_hits.out \
-r ejaculatory_sac_accepted_hits.out \
-r testes_accepted_hits.out \
-o Lycoriella_ingenua_scaffolds.rascaf
```

Then, to close gaps (Ns), I polished with Illumina reads once:
```
## prep Illumina reads for racon
cat *1.trimmed.fq.gz *2.trimmed.fq.gz > reads.combined.fastq.gz
zcat reads.combined.fastq.gz | tr -d " .*" | seqtk rename - reads. | awk '{print (NR%4 == 1 && NR%8 != 1) ? $0 ".2" : $0}' | awk '{print (NR%8 == 1) ? $0 ".1" : $0}' | gzip -c > reads.renamed.fastq.gz
# map and polish
minimap2 -t 60 -x sr ${genomename}.fa reads.renamed.fastq.gz > illum_vs_${genomename}.1.ctg.paf && racon -t 60 --include-unpolished --no-trimming reads.renamed.fastq.gz illum_vs_${genomename}.1.ctg.paf ${genomename}.fa > ${genomename}.i1.fa
```

## Gene annotations

Predicted gene sets were already publicly available for _Bradysia coprophila_, _Contarinia nasturtii_ and _Bradysia odoriphaga_. For all other species, I generated annotations using BRAKER2. Prior to this, I repeat-masked all genomes with Red:

```
Red -gnm ${base}/ -msk ${base}.red
```

Proportions of sequences masked are as follows:

| Species | Genome length | % masked |
| - | - | - |
| Aphidoletes aphidimyza | 191431736 | 33.5611% |
| Bradysia confinis |	338817951 | 31.213% |
| Bradysia desolata |	275897695 | 34.4701% |
| Bradysia pectoralis | 255331210 | 35.0162% |
| Conarinia rumicis	162540589 | 33.256% |
| Lycoriella agraria | 213236371 | 30.9621% |
| Lycoriella ingenua | 230556326 | 33.5557% |
| Mayetiola destructor | 185827756 | 43.6634% |
| Mayetiola hordei | 140973543 | 54.4051% |
| Phytosciara flavipes | 441741998 | 38.6771% |
| Sitodiplosis mosellana | 180693642 | 29.1561% |


To annotate, I used the OrthoDV v10 Diptera protein set as homology-based evidence. I annotated using the following command:

```
braker.pl \
--GENEMARK_PATH=~/software/gmes_linux_64 \
--AUGUSTUS_CONFIG_PATH=~/software/Augustus/config \
--AUGUSTUS_SCRIPTS_PATH=~/software/Augustus/scripts \
--PROTHINT_PATH=~/software/ProtHint/bin \
--genome=$GENOME \
--species=$GENOME.genome \
--gff3 \
--softmasking \
--verbosity=3 \
--cores 40 \
--useexisting \
--prot_seq=odb10v1_diptera_proteins.fasta \
--workingdir=$GENOME.braker
```




