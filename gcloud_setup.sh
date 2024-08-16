#!/bin/bash
set -e

mkdir -p ~/.local/bin
mkdir -p ~/Downloads
cd ~/Downloads || exit

# unzip
sudo apt install -y unzip

# bzip2
sudo apt install -y bzip2

#  SRA Toolkit https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar xvf sratoolkit.current-ubuntu64.tar.gz
mv sratoolkit.3.1.1-ubuntu64 ~/.local
cd ~/.local/bin || exit
ln -s ~/.local/sratoolkit.3.1.1-ubuntu64/bin/* .
cd || exit
mkdir -p ~/.cache/sra-toolkit

# seqkit
sudo apt update
sudo apt install -y seqkit

# pbzip2
sudo apt install -y pbzip2

# NCBI datasets
cd ~/Downloads || exit
wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets
mv datasets ~/.local/bin
chmod a+x ~/.local/bin/datasets

# samtools
sudo apt install -y samtools

# bowtie2
sudo apt install -y bowtie2

# bcftools
sudo apt install -y bcftools

