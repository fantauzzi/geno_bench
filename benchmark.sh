#!/bin/bash

set -e

get_timestamp() {
  date +%s
}

THREADS=16

# Accession number for the run that produced the reads
# ACC=ERR9466181 # Drosophila melanogaster
# ACC=SRR30202075 # Saccharomyces cerevisiae
# ACC=SRR17858636 # Encephalitozoon cuniculi
# ACC=SRR12007531 # Utricularia gibba
ACC=SRR5512131 #  Shigella flexneri 2a str. 301

# Accession number for the reference genome (NCBI datatbase)
# REF_ACC=GCF_000001215.4 # Drosophila melanogaster
# REF_ACC=GCF_000146045.2  # Saccharomyces cerevisiae
# REF_ACC=GCF_000091225.2 # Encephalitozoon cuniculi
# REF_ACC=GCA_002189035.1 # Utricularia gibba
REF_ACC=GCF_000006925.2 #  Shigella flexneri 2a str. 301

PLOIDY=1 # DON'T FORGET to set this correctly!

start_time=$(get_timestamp)
date

echo
echo "Experiment $ACC"
echo "Reference genome $REF_ACC"
echo "Ploidy is set to $PLOIDY"
echo
echo "Downloading reads"
# prefetch is way faster than just fastq-dump or fasterq-dump directly from remote
prefetch -v $ACC
# Convert the downloaded files to FASTQ
echo
echo "Converting downloaded reads to FASTQ format"
fasterq-dump -v --split-files $ACC

# Print information about the FASTQ files; must be done before compressing them as .bz2, cause
# seqkit doesn't read that compressed format
echo
# seqkit stats "$ACC"*.fastq.gz
seqkit stats -j $THREADS "$ACC"*.fastq

# Compress the FASTQ files, in parallel
echo
echo "Compressing FASTQ files"
# find . -maxdepth 1 -name "$ACC*.fastq" -print0 | xargs -0 -P "$THREADS" -I {} gzip {}
for file in ${ACC}*.fastq; do
  pbzip2 -p$THREADS -m2000 -v "$file"
done

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
READS1=$(ls "$ACC"*_1.fastq.bz2 | head -n 1)
READS2=$(ls "$ACC"*_2.fastq.bz2 | head -n 1)

# Index the reference genome to use with IGV
samtools faidx $REF_FILE

echo
echo "Indexing the reference for the aligner"
date
bowtie2-build --threads $THREADS -q $REF_FILE $REF_FILE

echo
echo "Running the aligner"
date
# bowtie2 --threads $THREADS --reorder -x $REF_FILE -1 $READS1 -2 $READS2 | samtools view --threads $THREADS -b > $ACC.bowtie2.bam
bowtie2 --threads $THREADS --reorder -x $REF_FILE -1 $READS1 -2 $READS2 | samtools sort --threads $THREADS -o $ACC.bowtie2.bam

echo
echo "Calling variants"
date
bcftools mpileup --threads $THREADS -O z -f $REF_FILE $ACC.bowtie2.bam | bcftools call --threads $THREADS --ploidy $PLOIDY -m -O z -o $ACC.bcftools.vcf

echo
date
end_time=$(get_timestamp)
elapsed_time=$((end_time - start_time))
echo "Time taken (elapsed): ${elapsed_time} sec."
