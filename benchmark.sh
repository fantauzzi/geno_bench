#!/bin/bash

set -e

get_timestamp() {
  date +%s
}

THREADS=$(nproc --all)

# Accession number for the run that produced the reads
# ACC=ERR9466181 # Drosophila melanogaster -longest time
# ACC=SRR29972059 # Saccharomyces cerevisiae -long time
# ACC=SRR17858636 # Encephalitozoon cuniculi -short time
ACC=SRR5512131 #  Shigella flexneri 2a str. 301 -shortest time

# Accession number for the reference genome (NCBI datatbase)
# REF_ACC=GCF_000001215.4 # Drosophila melanogaster
# REF_ACC=GCF_000146045.2  # Saccharomyces cerevisiae
# REF_ACC=GCF_000091225.2 # Encephalitozoon cuniculi
REF_ACC=GCF_000006925.2 #  Shigella flexneri 2a str. 301

PLOIDY=1 # DON'T FORGET to set this correctly!

start_time=$(get_timestamp)
date

echo
echo "Host name: $(hostname)"
echo "Experiment: $ACC"
echo "Reference genome: $REF_ACC"
echo "Ploidy is set to $PLOIDY"
echo "Number of threads: $THREADS"
echo
echo "Downloading reads"
# prefetch is way faster than just fastq-dump or fasterq-dump directly from remote
prefetch -v $ACC
# Convert the downloaded files to FASTQ
echo
echo "Converting downloaded reads to FASTQ format"
fasterq-dump -v --split-files $ACC

# Compress the FASTQ files, in parallel
echo
echo "Compressing FASTQ files"
# find . -maxdepth 1 -name "$ACC*.fastq" -print0 | xargs -0 -P "$THREADS" -I {} gzip {}
for file in ${ACC}*.fastq; do
  # pbzip2 -p$THREADS -m2000 -v "$file"
  bgzip -@ $THREADS $file
done

# Print information about the FASTQ files
echo
# seqkit stats "$ACC"*.fastq.gz
seqkit stats -j $THREADS "$ACC"*.fastq.gz

# Donwload the reference genome
rm -rf unzipped
rm -rf ncbi_dataset.zip
echo
echo "Obtaining the reference genome (if not already here)"
datasets download genome accession $REF_ACC
unzip ncbi_dataset.zip -d unzipped
cp "unzipped/ncbi_dataset/data/$REF_ACC/$REF_ACC"*.fna .
rm -rf unzipped

REF_FILE=$(ls "$REF_ACC"*.fna | head -n 1)
# READS1=$(ls "$ACC"*_1.fastq.bz2 | head -n 1)
READS1=$(ls "$ACC"*_1.fastq.gz | head -n 1)
# READS2=$(ls "$ACC"*_2.fastq.bz2 | head -n 1)
READS2=$(ls "$ACC"*_2.fastq.gz | head -n 1)

echo
date
echo "Compressing the reference genome"
# pbzip2 -p$THREADS -m2000 -v "$REF_FILE"
bgzip -@ $THREADS $REF_FILE
REF_FILE="$REF_FILE".gz

# Index the reference genome to use with IGV
samtools faidx $REF_FILE

echo
date
echo "Indexing the reference for the aligner"
bowtie2-build --threads $THREADS -q $REF_FILE $REF_FILE

echo
date
echo "Running the aligner"
bowtie2 --threads $THREADS -x $REF_FILE -1 $READS1 -2 $READS2 | samtools sort --threads $THREADS -o $ACC.bowtie2.bam

echo
date
echo "Computing stats for the alignment"
samtools flagstat -@ $THREADS $ACC.bowtie2.bam

echo
date
echo "Calling variants"
bcftools mpileup --threads $THREADS -O z -f $REF_FILE $ACC.bowtie2.bam | bcftools call --threads $THREADS --ploidy $PLOIDY -m -O z -o $ACC.bcftools.vcf

echo
date
end_time=$(get_timestamp)
elapsed_time=$((end_time - start_time))
echo
echo "Time taken (elapsed): ${elapsed_time} sec."

# TODO now do it with GATK
