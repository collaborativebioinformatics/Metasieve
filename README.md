# Metasieve
<img src="metasieve_logo.png" width="200">
## Overview of Metasieve

Mining environmental metagenomes for deeply novel biological sequence and structure

This is a metagenomic bioprospecting workflow designed to identify
environmental DNA and protein sequences that remain unexplained after
taxonomic, functional, and structural database searches.

## Team 
Dr. Todd Treangen, Rice University\
Dr. Jennifer Lu, The Johns Hopkins University\
Hiba Ben Aribi, Tunis El Manar University\
Felix Quintana, Rice University\
Natalie Kokroko, Rice University\
Jingyue Wu, Baylor College of Medicine\
Mohamed Abdelrahim


## Software Requirements
1. [MEGAHIT](https://github.com/voutcn/megahit)
2. [Kraken2](https://github.com/DerrickWood/kraken2/)
3. [Seqscreen](https://gitlab.com/treangenlab/seqscreen)
4. [ESMFold](https://esmatlas.com/)
5. [Foldseek](https://github.com/steineggerlab/foldseek)

## Workflow Overview


<img src="metasieve.png" width="450">

## Required Data
1. [Kraken2 Database: k2\_pluspf\_20260626.tar.gz](https://benlangmead.github.io/aws-indexes/k2): 
2. [Wastewater Metagenomic dataset 1](https://www.ncbi.nlm.nih.gov/sra/SRX34378765[accn]): 456.3G bases
3. [Wastewater Metagenomic dataset 2](https://www.ncbi.nlm.nih.gov/sra/SRX34378760[accn]): 488.7G bases
4. [Cacao Soil Metagenomics](https://www.ncbi.nlm.nih.gov/sra/SRX34351257[accn]): 25.6G bases

## Command Line
```
 prefetch \
    SRR000001 \
    --max-size 40g 
 fasterq-dump \
    --temp $TMPDIR \
    --split-files SRR000001
 megahit \
    -1 SRR000001_1.fastq.gz \
    -2 SRR000001_2.fastq.gz \
    -t $THREADS \
    -o SRR000001_megahit 
 kraken2 \
    --db $KRAKEN_DB \
    --threads 16 \
    --confidence 0.01 \
    --report SRR000001_scaffolds.k2report \
    --unclassified-out SRR000001_k2unclassified.fa \
    --output SRR000001_scaffolds.k2 
    SRR000001_scaffolds.fasta
 seqscreen \
    --fasta SRR000001_k2unclassified.fa \
    --databases SeqScreenDB_23.4 \
    --working seqscreen_out \
    --threads 16 \
```

## Ongoing/Upcoming Projects (Updated 9AM CST 2026/08/27)
1. Run full pipeline (start to finish) on cacao soil test-set
2. Run full pipeline (start to finish) on Wastewater dataset
3. Incorporate human genome removal (bowtie2 against T2T) pre-assembly/classification
4. Compare assembly pre/post Kraken2 classification
5. Compare megahit vs. metaspades vs. ggcat for assembly 
