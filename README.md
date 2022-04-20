# Yleaf: software for human Y-chromosomal haplogroup inference from next generation sequencing data

#### Arwin Ralf, Diego Montiel Gonzalez, Kaiyin Zhong and Manfred Kayser

### Department of Genetic Identification 
#### Erasmus MC University Medical Centre Rotterdam, The Netherlands

## Requirements

    Operating system: Linux only. Tested on Ubuntu 16.04LTS, but should also work on newer version of Ubuntu. It should be easy to made it work on other Linux distributions. 
    Python3 with pandas >= 0.23.0 and numpy >= 1.14.3
    Internet connection during installation (for downloading and extracting hg19/hg38 reference genome).
    Data storage: For installation we recommend a storage capacity of > 10 GB. 

## Installation

You can either build a conda environment with all dependencies

```bash
# create the conda environment
conda create --name yleaf_env python=3.6
# activat the environment
conda activate yleaf_env
# install required python libraries
conda install pandas
conda install numpy
# install Burrows-Wheeler Aligner for fastQ files
conda install -c bioconda bwa
# install SAMtools
conda install -c bioconda samtools
# clone the yleaf repository
git clone https://github.com/genid/Yleaf.git

# finally if you want to use FASTQ files use the provided install.py
script to download required genome data
cd Yleaf
python install.py
```      
or manually install everything
```bash
# install python and libraries
apt-get install python3.6
pip3 install pandas
pip3 install numpy
# install Burrows-Wheeler Aligner for FASTQ files
sudo apt-get install bwa
# install SAMtools
wget https://github.com/samtools/samtools/releases/download/1.4.1/
samtools-1.4.1.tar.bz2 -O samtools.tar.bz2
tar -xjvf samtools.tar.bz2 3. cd samtools-1.4.1/
./configure 5. make
make install
# clone the yleaf repository
git clone https://github.com/genid/Yleaf.git

# finally if you want to use FASTQ files use the provided install.py
script to download required genome data
cd Yleaf
python install.py
```
## Usage and examples

In version 3.0 we switched to using YFull (v10.01) for the underlying tree structure of the haplogroups.
 This also means that predictions are a bit different compared to earlier versions.
### Yleaf: FASTQ (raw reads)
    
    python Yleaf.py -fastq raw_reads.fastq -f reference_indexed/genome.fasta -pos Position_files/new_position_files/[hg19.txt/hg38.txt] -out out -r 1 -q 20 -b 90 -t 4
        
### Yleaf: BAM or CRAM format
    
    python Yleaf.py -bam file.bam -pos Position_files/new_position_files/[hg19.txt/hg38.txt] -out out -r 1 -q 20 -b 90 

    python Yleaf.py -cram file.cram -pos Position_files/new_position_files/[hg19.txt/hg38.txt] -out out -r 1 -q 20 -b 90  -f genome.fasta

### Haplogroup prediction (w/output generated from Yleaf)

    python predict_haplogroup.py -input Output_files/ -out output.hg

## Additional information

For a more comprehensive manual please have a look at the [yleaf_manual](yleaf_manual.pdf).

If you have a bug to report or a question about installation consider sending an email to 
 b.vanwersch at erasmusmc.nl or create an issue on GitHub.


### References and Supporting Information
A. Ralf, et al., Yleaf: software for human Y-chromosomal haplogroup inference from next generation sequencing data (2018).

https://academic.oup.com/mbe/article/35/7/1820/4993044

