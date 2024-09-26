configfile: 'config.yaml'
locals().update(config)
REF_FILENAME ='refs/' + REF_ACC+'.zip'

rule all:
    input:
        f'{REF_DATA_DIR}/{REF_ACC}/{REF_ACC}.fna.gz.fai',
        f'{REF_DATA_DIR}/{REF_ACC}/{REF_ACC}.fna.gz.gzi',
        f'{BAM_DIR}/{ACC}-{REF_ACC}.bowtie2.bam'

rule clean:
    shell:
        'rm -rf bam fastq refs sra'

rule download_reads:
    output:
        sra = f'{SRA_DIR}/{ACC}/{ACC}.sra'
    shell:
        "prefetch -v {ACC} -O {SRA_DIR}"

rule convert_sra_to_fastq:
    input:
        sra = rules.download_reads.output.sra
    output:
        fastq1 = f'{FASTQ_DIR}/{ACC}_1.fastq',
        fastq2 = f'{FASTQ_DIR}/{ACC}_2.fastq'
    shell:
        "fasterq-dump -v -O {FASTQ_DIR} --split-files {input.sra}"

rule compress_fastq1_file:
    output:
        fastq1_gz = f'{FASTQ_DIR}/{ACC}_1.fastq.gz'
    input:
        rules.convert_sra_to_fastq.output.fastq1
    shell:
        "bgzip -@ {THREADS} {input};"
        "seqkit stats -j {THREADS} {output}" # Send this to a log file, instead of cluttering the output of snakemake

rule compress_fastq2_file:
    output:
        fastq2_gz = f'{FASTQ_DIR}/{ACC}_2.fastq.gz'
    input:
        rules.convert_sra_to_fastq.output.fastq2
    shell:
        "bgzip -@ {THREADS} {input};"
        "seqkit stats -j {THREADS} {output}" # Send this to a log file, instead of cluttering the output of snakemake

rule download_ref_genome:
    output:
        ref = f'{REF_FILENAME}'
    shell:
        "datasets download genome accession {REF_ACC} --filename {output.ref}"

rule unzip_ref_genome:
    output:
        fna_file = f'{REF_DATA_DIR}/{REF_ACC}/{REF_ACC}.fna'
    input:
        rules.download_ref_genome.output.ref
    shell:
        "unzip -u {input} -d {REFS_DIR}; "
        'mv "$(ls {REF_DATA_DIR}/{REF_ACC}/{REF_ACC}*.fna)" {output.fna_file}'

rule bgzip_fna:
    input:
        rules.unzip_ref_genome.output
    output:
        bgzipped = rules.unzip_ref_genome.output.fna_file + '.gz'
    shell:
        "bgzip -f -@ {THREADS} {input}"

rule index_reference_for_igv:
    input:
        rules.bgzip_fna.output
    output:
        rules.bgzip_fna.output.bgzipped + '.fai',
        rules.bgzip_fna.output.bgzipped + '.gzi'
    shell:
        'samtools faidx {rules.bgzip_fna.output}'

rule index_reference_for_bowtie2:
    input:
        rules.bgzip_fna.output
    output:
        rules.bgzip_fna.output.bgzipped + '.1.bt2',  # TODO use expand() here
        rules.bgzip_fna.output.bgzipped + '.2.bt2',
        rules.bgzip_fna.output.bgzipped + '.3.bt2',
        rules.bgzip_fna.output.bgzipped + '.4.bt2',
        rules.bgzip_fna.output.bgzipped + '.rev.1.bt2',
        rules.bgzip_fna.output.bgzipped + '.rev.2.bt2'
    shell:
        'bowtie2-build --threads {THREADS} -q {rules.bgzip_fna.output.bgzipped} {rules.bgzip_fna.output.bgzipped}'

rule align_reads:
    input:
        rules.bgzip_fna.output,
        rules.index_reference_for_bowtie2.output,
        rules.compress_fastq1_file.output,
        rules.compress_fastq2_file.output
    output:
        bam = f'{BAM_DIR}/{ACC}-{REF_ACC}.bowtie2.bam'
    shell:
        'bowtie2 --threads {THREADS} -x {rules.bgzip_fna.output} -1 {rules.compress_fastq1_file.output} -2 {rules.compress_fastq2_file.output} | samtools sort --threads {THREADS} -o {output.bam}'
