
SAMPLES = ["CR_H3K27ac_cranial_1", "CR_H3K27ac_cranial_2", "CR_H3K27ac_trunk_1", "CR_H3K27ac_trunk_2", "CR_Smad23_1", "CR_Smad23_2", "CR_bCat_1", "CR_bCat_2", "CR_IgG"]

rule all:
	input:
		expand("fastq/{sample}_R1_fastqc.html", sample=SAMPLES),
		expand("trimmed_fastq/trimmed_{sample}_R1.fastq.gz", sample=SAMPLES),
		expand("trimmed_fastq/trimmed_{sample}_R2.fastq.gz", sample=SAMPLES),
		expand("logs/cutadapt/{sample}.log", sample=SAMPLES),
		expand("bam/galGal6/{sample}_galGal6.bam", sample=SAMPLES),
		expand("bam/galGal6/{sample}_galGal6.bam.bai", sample=SAMPLES),
		expand("peaks/{sample}_peaks.narrowPeak", sample=SAMPLES),
		expand("bw/{sample}_galGal6.bw", sample=SAMPLES)


rule fastqc:
	input:
		R1="fastq/{sample}_R1.fastq.gz",
		R2="fastq/{sample}_R2.fastq.gz"
	output:
		R1="fastq/{sample}_R1_fastqc.html",
		R2="fastq/{sample}_R2_fastqc.html"
	threads: 8
	shell:
		"fastqc {input} -t {threads}"


rule cutadapt:
	input:
		R1="fastq/{sample}_R1.fastq.gz",
		R2="fastq/{sample}_R2.fastq.gz"
	output:
		R1="trimmed_fastq/trimmed_{sample}_R1.fastq.gz",
		R2="trimmed_fastq/trimmed_{sample}_R2.fastq.gz"
	log:
		"logs/cutadapt/{sample}.log"
	threads: 8
	shell:
		"cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
		-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
		--minimum-length=25 \
		--cores={threads} \
		{input.R1} \
		{input.R2} \
		-o {output.R1} \
		-p {output.R2} 2> {log}"


rule bowtie2:
	input:
		R1="trimmed_fastq/trimmed_{sample}_R1.fastq.gz",
		R2="trimmed_fastq/trimmed_{sample}_R2.fastq.gz"
	output:
		bam="bam/galGal6/{sample}_galGal6.bam"
	log:
		"logs/bowtie2/{sample}_galGal6.log"
	threads: 16
	shell:
		"bowtie2 --local --very-sensitive-local \
		--no-unal --no-mixed --no-discordant \
		-x /data/Megan/genome_data/galGal6/galGal6_bt2/galGal6_bt2 \
		-I 10 -X 1000 \
		-1 {input.R1} \
		-2 {input.R2} \
		--threads {threads} \
		2> {log} | samtools view -F 780 -f 2 -bh - | samtools sort -T {wildcards.sample}_galGal6 - -o {output.bam}"


rule index:
	input:
		"bam/galGal6/{sample}_galGal6.bam"
	output:
		"bam/galGal6/{sample}_galGal6.bam.bai"
	shell:
		"samtools index {input}"


rule callpeak:
	input:
		"bam/galGal6/{sample}_galGal6.bam"
	output:
		"peaks/{sample}_peaks.narrowPeak"
	log:
		"logs/macs2/{sample}.log"
	shell:
		"macs2 callpeak -t {input} -n {wildcards.sample} -f BAMPE -g 1e9 -q 0.05 --call-summits --outdir ./peaks 2> {log}"


rule bamCoverage:
	input:
		bam="bam/galGal6/{sample}_galGal6.bam", 
		bai="bam/galGal6/{sample}_galGal6.bam.bai"
	output:
		"bw/{sample}_galGal6.bw"
	threads: 8
	shell:
		"bamCoverage --bam {input.bam} \
		--outFileName {output} \
		--outFileFormat bigwig \
		--binSize 5 \
		--numberOfProcessors {threads} \
		--normalizeUsing RPKM \
		--extendReads"



