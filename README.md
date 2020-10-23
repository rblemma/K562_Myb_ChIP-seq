# K562_Myb_ChIP-seq
Snakemake workflow for processing in-house Myb ChIP-seq data (GEO Accession: GSE124541)

## Raw reads 
* Raw reads can be obtained from GSE124541 and put in the 'data/raw_reads/' folder
* PhiX genome are placed in 'data/PhiX_fasta' folder
* The Snakefile used to process the c-Myb ChIP-seq data can be found in 'results/' folder

## Softwares
* snakemake
* trim-galore
* bbmap
* bwa
* samtools
* MACS2
* intervene

### Additional softwares used in the paper but not in the Snakefile
* deeptools
* HOMER
* bedtools

