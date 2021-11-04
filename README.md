# NT-seq
This script is designed for data analysis of NT-seq.

## Requirement: 
Python3 \
Bowtie2 \
Cutadapt \
Snakemake \
Samtools \
Igvtools \
Biopython

## Usage: 
Create conda environment:\
`conda create -f NT-seq.yml`\
Activate conda environment:\
`conda activate NT-seq`\
Run workflow with snakemake:\
`snakemake -s NT-seq-workflow_PE.py -j 10`

## Expected output: 
CSV files containing base count at each position of reference genome for each sample, which can be used to calculate A to G ratio and C to T ratio for 6mA, 4mC, and 5mC motif calling.

## Test data: 
random sampled H. pylori NT-seq data.

## Expected run time for demo
About ten minutes using ten threads of Intel(R) Xeon(R) Gold 5120 CPU @ 2.20GHz

## Notes: 
To apply the script to other data, please modify the ref_genome and sample name in snakemake scripts.
