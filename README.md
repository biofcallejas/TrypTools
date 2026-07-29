# TrypTools
These Python scripts identify, confirm, and extract SL sequences from RNA-seq reads. 
The output also includes trimmed (SL) fastq reads.

## Installation 

> [!NOTE]
Although no installation is needed for the main script, dependencies/python packages are needed.

The script has been tested on Python 3.8.12, Mac (Intel) and Linux.

Python package needed (double-check it's installed and available): ```ProcessPoolExecutor```

> [!IMPORTANT]
> Before running the scripts, the following dependencies need to be installed and available on the path.

```
Bowtie2 (tested on 2.4.1)
gffread (tested on v0.12.7)
```

## Overview

The script is divided into three main steps:

```
1.- SL sequence is searched (3 and 5 ends )
2.- RNAseq reads are mapped to a reference genome
3.- SL sequences are confirmed not to be part of the same genomic region as the rest of the read (true positives).
```

## Mandatory arguments:

```
-s,  --spliced,  The reference SL sequence in fasta format
-g,  --genome,  Genome reference, fasta format
-i,  --gindex,  Prefix for the Bowtie genome index

for paired-end reads:
-r1,  --fastq1,  Fastq reads mate1; fragment r1
-r2,  --fastq2,  Fastq reads mate2; fragment r2

For single-end reads:
-q,  --fastq,  Fastq reads, single end
```
## Optional arguments:

```
-t,  --threads, Number of threads to use, default=1
-p,  --sprefix, Use this parameter if you're running multiple jobs on the same path, the script creates multiple temporary files that will cause redundant files, default=search_sl_trimmer
-m,  --mismatch, Maximum number of mismatches in the SL, default=6
-l,  --slength, Minimum length of the SL, default=8
```

## Run example:

**For paired-end** 

```
python search_SL_TC_paired-end_v2.py -s trypanosoma.fa -r1 SRR25010101_1.fastq -r2 SRR25010101_2.fastq -g TcruziBerenice.fasta -i TcruziBerenice_idx -t 20 -p SRR25010101
```

**For single-end** 

```
python search_SL_TC_single-end_v5.py -s trypanosoma.fa -q SRR28628196.fastq -g TcruziBerenice.fasta -i TcruziBerenice_idx -t 20 -p SRR28628196
```

## Output example:

SL search in fastq file has started...
Time elapsed: 5.97 minutes
Potential Leader hits: 1102212

Bowtie2 alignment has started has started...

1102212 reads; of these:
  1102212 (100.00%) were paired; of these:
    314112 (28.50%) aligned concordantly 0 times
    547554 (49.68%) aligned concordantly exactly 1 time
    240546 (21.82%) aligned concordantly >1 times
    ----
    314112 pairs aligned concordantly 0 times; of these:
      66801 (21.27%) aligned discordantly 1 time
    ----
    247311 pairs aligned 0 times concordantly or discordantly; of these:
      494622 mates make up the pairs; of these:
        332369 (67.20%) aligned 0 times
        59464 (12.02%) aligned exactly 1 time
        102789 (20.78%) aligned >1 times
84.92% overall alignment rate

Bowtie2 alignment has finished...

SAM analisis has started...
SAM analisis has finished...

Valid leader hits: 925575

