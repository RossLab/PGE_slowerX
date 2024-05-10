# 1. Data

Genomes for the following species were downloaded from the FTP directories with `wget`. There are two BioProjects for _B. coprophila_ because one is for the X' chromosome and the other is for the core genome.

| species | BioProject |
| - | - |
| Aphidoletes aphidimyza | PRJNA634691 |
| Bradysia coprophila | PRJNA953429, PRJNA291918 |
| Bradysia odoriphaga | PRJNA612767 |
| Contarinia nasturtii | PRJNA612767 |
| Dilophus febrilis | PRJEB64321 |
| Mayetiola destructor | PRJNA45867 |
| Phytosciara flavipes | PRJNA369255 |
| Sitodiplosis mosellana | PRJNA720212 |

The following WGS libraries were downloaded from NBCI using the SRA toolkit, i.e.:

```
cat SRA_list.txt | while read SRA_ID
do
	~/software/sra_tools/sratoolkit.3.0.7-ubuntu64/bin/prefetch --max-size 200g $SRA_ID
	~/software/sra_tools/sratoolkit.3.0.7-ubuntu64/bin/fasterq-dump $SRA_ID
done 
```

| species | SRA |
| - | - |
| Aphidoletes aphidimya | SRR14087947 |
| Bradysia odoriphaga | SRR13756224 |
| Contarinia nasturtii | SRR11470102 |
| Contarinia rumicis | SRR10116830 |
| Dilophus febrilis | ERR12148423 |
| Dilophus femoratus | ERR12148424 |
| Mayetiola destructor | SRR1738190 |
| Mayetiola hordei | SRR136493 |
| Sitodiplosis mosellana | SRR14251190 |


We also generated the following data:

- 1 _Bradysia confinis_ male WGS library
- 1 _Bradysia pectoralis_ male WGS library
- 1 _Bradysia desolata_ male WGS library
- 1 _Lycoriella ingenua_ male WGS library
- 1  _Lycoriella agraria_ female WGS library
- 11 _Bradysia coprophila_ female WGS libraries
- 5 _Lycoriella ingenua_ female WGS libraries
- 18 _Bradysia coprophila_ RNAseq libraries
- 18 _Lycoriella ingenua_ RNAseq libraries

All libraries were trimmed with fastp and quality-checked with FastQC:

```
for file in $(ls *_1.fastq.gz)
do
	base=$(basename $file "_1.fastq.gz")
  	fastqc ${base}_1.fastq.gz ${base}_2.fastq.gz
  	fastp -i ${base}_1.fastq.gz -I ${base}_2.fastq.gz -o ${base}_1.trimmed.fastq.gz -O ${base}_2.trimmed.fastq.gz && \
  	rm ${base}_1.fastq.gz ${base}_2.fastq.gz
done
```




















