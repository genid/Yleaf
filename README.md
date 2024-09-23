# Yleaf: software for human Y-chromosomal haplogroup inference from next generation sequencing data

#### Arwin Ralf, Diego Montiel Gonzalez, Kaiyin Zhong and Manfred Kayser

### Department of Genetic Identification 
#### Erasmus MC University Medical Centre Rotterdam, The Netherlands

## Requirements

    Operating system: Linux only. 
    Internet connection: when running for the first time for downloading the reference genome. Alternatively you 
                         can configure your own references.
    Data storage: For installation we recommend a storage capacity of > 8 GB. 

## Installation

The easiest way to get Yleaf up and running is by using a conda environment. 

```bash
# first clone this repository to get the environment_yleaf.yaml
git clone https://github.com/genid/Yleaf.git
cd Yleaf
# create the conda environment from the .yaml the environment will be called yleaf
conda env create --file environment_yleaf.yaml
# activate the environment
conda activate yleaf
# pip install the cloned yleaf into your environment. Using the -e flag allows you to modify the config file in your cloned folder
pip install -e .

# verify that Yleaf is installed correctly. You can call this command from any directory on your system
Yleaf -h 
```      
or manually install everything
```bash
# install python and libraries
apt-get install python3.6
pip3 install pandas
pip3 install numpy
# install Burrows-Wheeler Aligner for FASTQ files
sudo apt-get install minimap2 
# install SAMtools
wget https://github.com/samtools/samtools/releases/download/1.4.1/
samtools-1.4.1.tar.bz2 -O samtools.tar.bz2
tar -xjvf samtools.tar.bz2
cd samtools-1.4.1/
./configure make
make install
# clone the yleaf repository
git clone https://github.com/genid/Yleaf.git
# pip install the yleaf repository
cd Yleaf
pip install -e .

# verify that Yleaf is installed correctly. You can call this command from any directory on your system
Yleaf -h 
```
After installation you can navigate to the yleaf/config.txt folder and add custom paths for the files listed there. This will make sure that Yleaf does not download the files on the first go or downloads the files in the provided location. This allows you to use a custom reference if you want. Please keep in mind that custom reference files might cause other issues or give problems in combination with already existing data files. Positions are based on either hg38 or hg19.

## Usage and examples
Here follow some minimal working examples of how to use Yleaf with different input files. There are additional options
that can be used to tune how strict Yleaf is as well as options to get private mutations as well as a graph showing 
the positioning of predicted haplogroups of all your samples in the Haplogroup tree.

_Note: In version 3.0 we switched to using YFull (v10.01) for the underlying tree structure of the haplogroups.
 This also means that predictions are a bit different compared to earlier versions._
### Yleaf: FASTQ (raw reads)

    Yleaf -fastq raw_reads.fastq -o fastq_output --reference_genome hg38
        
### Yleaf: BAM or CRAM format
    Yleaf -bam file.bam -o bam_output --reference_genome hg19 
    Yleaf -cram file.bam -o cram_output --reference_genome hg38 

### With drawing predicted haplogroups in a tree and showing all private mutations

    Yleaf -bam file.bam -o bam_output --reference_genome hg19 -dh -p

## Additional information

For a more comprehensive manual please have a look at the [yleaf_manual](yleaf_manual.pdf).

If you have a bug to report or a question about installation consider sending an email to 
 a.ralf at erasmusmc.nl or create an issue on GitHub.


### References and Supporting Information
A. Ralf, et al., Yleaf: software for human Y-chromosomal haplogroup inference from next generation sequencing data (2018).

https://academic.oup.com/mbe/article/35/5/1291/4922696

