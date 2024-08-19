ACC='SRR5512131'
REF_ACC='GCF_000006925.2' #  Shigella flexneri 2a str. 301
SRA_DIR = 'sra'
FASTQ_DIR = 'fastq'
REF_FILENAME ='refs/' + REF_ACC+'.zip'
REFS_DIR = 'refs'
REF_DATA_DIR = 'refs/ncbi_dataset/data'
BAM_DIR = 'bam'
THREADS = 16

rule all:
    input:
        igv_index_fai = expand('{ref_data_dir}/{ref_acc}/{ref_acc}.fna.gz.fai', ref_data_dir=REF_DATA_DIR, ref_acc=REF_ACC),
        igv_index_gzi = expand('{ref_data_dir}/{ref_acc}/{ref_acc}.fna.gz.gzi', ref_data_dir=REF_DATA_DIR, ref_acc=REF_ACC),
        bam_file = expand("{bam_dir}/{acc}-{ref_acc}.bowtie2.bam",acc=ACC, ref_acc=REF_ACC, bam_dir=BAM_DIR)

rule download_reads:
    output:
        sra = expand("{sra_dir}/{acc}/{acc}.sra", sra_dir=SRA_DIR, acc=ACC)
    shell:
        "prefetch -v {ACC} -O {SRA_DIR}"

rule convert_sra_to_fastq:
    output:
        fastq1 = expand("{fastq_dir}/{acc}_1.fastq", acc=ACC, fastq_dir=FASTQ_DIR),
        fastq2 = expand("{fastq_dir}/{acc}_2.fastq", acc=ACC, fastq_dir=FASTQ_DIR)
    input:
        rules.download_reads.output
    shell:
        "fasterq-dump -v -O {FASTQ_DIR} --split-files {input}"

rule compress_fastq1_file:
    output:
        fastq1_gz = expand("{fastq_dir}/{acc}_1.fastq.gz", acc=ACC, fastq_dir=FASTQ_DIR),
    input:
        rules.convert_sra_to_fastq.output.fastq1
    shell:
        "bgzip -@ {THREADS} {input};"
        "seqkit stats -j {THREADS} {output}" # Send this to a log file, instead of cluttering the output of snakemake

rule compress_fastq2_file:
    output:
        fastq2_gz = expand("{fastq_dir}/{acc}_2.fastq.gz", acc=ACC, fastq_dir=FASTQ_DIR),
    input:
        rules.convert_sra_to_fastq.output.fastq2
    shell:
        "bgzip -@ {THREADS} {input};"
        "seqkit stats -j {THREADS} {output}" # Send this to a log file, instead of cluttering the output of snakemake

rule download_ref_genome:
    output:
        ref = expand("{ref_filename}", ref_filename=REF_FILENAME)
    shell:
        "datasets download genome accession {REF_ACC} --filename {output.ref}"

rule unzip_ref_genome:
    output:
        fna_file = expand('{ref_data_dir}/{ref_acc}/{ref_acc}.fna', ref_data_dir=REF_DATA_DIR, ref_acc=REF_ACC)
    input:
        rules.download_ref_genome.output
    shell:
        "unzip -u {input} -d {REFS_DIR}; "
        'mv "$(ls {REF_DATA_DIR}/{REF_ACC}/{REF_ACC}*.fna)" {REF_DATA_DIR}/{REF_ACC}/{REF_ACC}.fna'


rule bgzip_fna:
    input:
        rules.unzip_ref_genome.output
    output:
        rules.unzip_ref_genome.output.fna_file[0] + '.gz'
    shell:
        "bgzip -@ {THREADS} {input}"

rule index_reference_for_igv:
    input:
        rules.bgzip_fna.output
    output:
        rules.bgzip_fna.output[0] + '.fai',
        rules.bgzip_fna.output[0] + '.gzi'
    shell:
        'samtools faidx {rules.bgzip_fna.output}'

rule index_reference_for_bowtie2:
    input:
        rules.bgzip_fna.output
    output:
        rules.bgzip_fna.output[0] + '.1.bt2',
        rules.bgzip_fna.output[0] + '.2.bt2',
        rules.bgzip_fna.output[0] + '.3.bt2',
        rules.bgzip_fna.output[0] + '.4.bt2',
        rules.bgzip_fna.output[0] + '.rev.1.bt2',
        rules.bgzip_fna.output[0] + '.rev.2.bt2'
    shell:
        'bowtie2-build --threads {THREADS} -q {rules.bgzip_fna.output[0]} {rules.bgzip_fna.output[0]}'

rule align_reads:
    input:
        rules.bgzip_fna.output,
        rules.index_reference_for_bowtie2.output,
        rules.compress_fastq1_file.output,
        rules.compress_fastq2_file.output
    output:
        expand("{bam_dir}/{acc}-{ref_acc}.bowtie2.bam",acc=ACC, ref_acc=REF_ACC, bam_dir=BAM_DIR)
    shell:
        'bowtie2 --threads {THREADS} -x {rules.bgzip_fna.output} -1 {rules.compress_fastq1_file.output} -2 {rules.compress_fastq2_file.output} | samtools sort --threads {THREADS} -o {BAM_DIR}/{ACC}-{REF_ACC}.bowtie2.bam'
