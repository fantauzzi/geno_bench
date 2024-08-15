module load  Bowtie2
module load  SAMtools/1.12-GCC-10.3.0


cd ~/Downloads || exit
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-ubuntu64.tar.gz
tar xvf sratoolkit.3.1.1-ubuntu64.tar.gz
mv mv sratoolkit.3.1.1-ubuntu64 ~/.local
cd ~/.local/bin || exit
ln -s ../sratoolkit.3.1.1-ubuntu64/bin/* .

cd ~/Downloads || exit
wget https://github.com/shenwei356/seqkit/releases/download/v2.8.2/seqkit_linux_amd64.tar.gz
tar xvf seqkit_linux_amd64.tar.gz
mv seqkit ~/.local/bin

wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets
mv datasets ~/.local/bin