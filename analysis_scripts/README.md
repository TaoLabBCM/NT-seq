This script is designed for data analysis of NT-seq.

Requirement:
Python3
Bowtie2
Cutadapt
Snakemake
Samtools
Igvtools
Biopython

Usage:
conda create -f NT-seq.yml
conda activate NT-seq
snakemake -s nitrite-seq-workflow_PE.py -j 10

Expected output:
CSV files containing base count at each position of reference genome for each sample.

Test data:
random sampled H. pylori NT-seq data.

Notes:
To apply the script to other data, please modify the ref_genome and sample name in snakemake scripts.