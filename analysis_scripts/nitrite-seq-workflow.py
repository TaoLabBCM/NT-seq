SAMPLES = ["Zymo-microbial-std-PCR", "Zymo-microbial-std-native"]
Read = ["R1", "R2"]
Strand = ["fwd", "rev"]
ruleorder: trimming>dedupe>fastq_AT_convert>bowtie2 > sam_to_original_pe > samtools_calmd > sam_tag_filtering > samtools_sort_index > igvtools
rule all:
    input:
        expand("{sample}_trimmed_{read}.fastq", sample=SAMPLES, read=Read),
        expand("{sample}_trimmed_dedupe_{read}.fastq", sample=SAMPLES, read=Read),
        expand("{sample}_trimmed_dedupe_{read}_AT_only.fastq", sample=SAMPLES, read=Read),
        expand("{sample}_AT_only.sam", sample=SAMPLES),
        expand("{sample}_original_reads_{strand}.sam", sample=SAMPLES, strand=Strand),
        expand("{sample}_original_reads_{strand}.sorted.bam", sample=SAMPLES, strand=Strand),
        expand("{sample}_original_reads_{strand}_baq.sam", sample=SAMPLES, strand=Strand),
        expand("{sample}_original_reads_{strand}_baq_filtered.sam",sample=SAMPLES, strand=Strand),
        expand("{sample}_original_reads_{strand}_baq_filtered.sorted.bam",sample=SAMPLES, strand=Strand),
        expand("readcount/{sample}_original_reads_{strand}_baq_filtered.wig",sample=SAMPLES, strand=Strand),


rule trimming:
    input:
        r1 = "{sample}_R1.fastq.gz",
        r2 = "{sample}_R2.fastq.gz"
    output:
        r1 = "{sample}_trimmed_R1.fastq",
        r2 = "{sample}_trimmed_R2.fastq"
        
    log:
        "logs/cutadpt/{sample}.log"
    params:
        TRIM_OPTS = "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -m 10"
    threads:
        32
    shell:
        """
        cutadapt {params.TRIM_OPTS} -j {threads} -o {output.r1} -p {output.r2} {input.r1} {input.r2}  1>{log}
        """
rule dedupe:
    input:
        r1 = "{sample}_trimmed_R1.fastq",
        r2 = "{sample}_trimmed_R2.fastq"
    output:
        r1 = "{sample}_trimmed_dedupe_R1.fastq",
        r2 = "{sample}_trimmed_dedupe_R2.fastq"
    log:
        "logs/dedupe/{sample}.log"
    shell:
        """
        clumpify.sh in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} dedupe 
        """
rule fastq_AT_convert:
    input:
        "{sample}_trimmed_dedupe_{read}.fastq"
    output:
        "{sample}_trimmed_dedupe_{read}_AT_only.fastq"
    shell:
        """
        python3 fastq_to_AT_only.py {input}
        """

rule bowtie2:
    input:
        r1 = "{sample}_trimmed_dedupe_R1_AT_only.fastq",
        r2 = "{sample}_trimmed_dedupe_R2_AT_only.fastq"
    params:
        index = "/storage/taowu/home/xuwenl/ref_genome/Zymo_microbial/bowtie2_index_AT_only_new/zymo_microbial_std"
    output:
        "{sample}_AT_only.sam"
    log:
        "logs/bowtie/{sample}.log"
    threads: 30
    shell:
        "bowtie2 -x {params.index} -p {threads} -1 {input.r1} -2 {input.r2} -S {output} 2>{log}"

rule sam_to_original_pe:
    input:
        sam = "{sample}_AT_only.sam",
        r1 = "{sample}_trimmed_dedupe_R1.fastq",
        r2 = "{sample}_trimmed_dedupe_R2.fastq"
    output:
        "{sample}_original_reads_fwd.sam",
        "{sample}_original_reads_rev.sam"
    shell:
        """
        python3 sam_to_original_pe.py {input.sam} {input.r1} {input.r2} {output[0]} {output[1]} 
        """
rule samtools_calmd:
    input:
        "{sample}_original_reads_{strand}.sam"
    output:
        "{sample}_original_reads_{strand}.sorted.bam",
        "{sample}_original_reads_{strand}_baq.sam"
    log:
        "logs/samtools_calmd/{sample}_{strand}.log"
    threads:
        10
    shell:
        """
        samtools sort {input} > {output[0]}
        samtools calmd -@ {threads} -Ar {output[0]} /storage/taowu/home/xuwenl/ref_genome/Zymo_microbial/zymo_microbial_std_new.fasta > {output[1]} 2>{log}
        """
rule sam_tag_filtering:
    input:
        fwd = "{sample}_original_reads_fwd_baq.sam",
        rev = "{sample}_original_reads_rev_baq.sam"
    output:
        fwd = "{sample}_original_reads_fwd_baq_filtered.sam",
        rev = "{sample}_original_reads_rev_baq_filtered.sam"

    shell:
        """
        python3 sam_tag_filtering.py {input.fwd} {output.fwd}
        python3 sam_tag_filtering.py {input.rev} {output.rev}
        """
rule samtools_sort_index:
    input:
        "{sample}_original_reads_{strand}_baq_filtered.sam"
    output:
        "{sample}_original_reads_{strand}_baq_filtered.sorted.bam"

    shell:
        """
        samtools sort -@ 20 {input} > {ouput[0]}
        samtools index {ouput[0]}
        """
rule igvtools:
    input:
        bam = "{sample}_original_reads_{strand}_baq_filtered.sorted.bam"
    output:
        "readcount/{sample}_original_reads_{strand}_baq_filtered.wig"
    shell:
        """
        igvtools count --bases -w 1 -e 0 {input.bam} {output} /storage/taowu/home/xuwenl/ref_genome/Zymo_microbial/zymo_microbial_std_new.fasta
        """
