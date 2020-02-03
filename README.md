# Yleaf: software for human Y-chromosomal haplogroup inference from next generation sequencing data

#### Arwin Ralf, Diego Montiel Gonzalez, Kaiyin Zhong and Manfred Kayser

### Department of Genetic Identification 
#### Erasmus MC University Medical Centre Rotterdam, The Netherlands

## Requirements

    Operating system: Linux only. Tested on Ubuntu 16.04LTS, but should also work on newer version of Ubuntu. It should be easy to made it work on other Linux distributions. 
    Python3 with pandas >= 0.23.0 and numpy >= 1.14.3
    Internet connection during installation (for downloading and extracting hg19 reference genome).
    Data storage: For installation we recommend a storage capacity of > 10 GB. 

## Installation

1. Install dependencies, you can skip this step if these packages are already installed on your system
            
            apt-get install python3.6            
            pip3 install pandas --upgrade
            pip3 install numpy  --upgrade            
            apt-get install bwa            
            

	SAMtools: We recommend the newests versions of SAMtools (e.g. > 1.4.1)

            1. wget https://github.com/samtools/samtools/releases/download/1.4.1/samtools-1.4.1.tar.bz2 -O samtools.tar.bz2
            2. tar -xjvf samtools.tar.bz2 
            3. cd samtools-1.4.1/
            4. ./configure
            5. make
            6. make install

## Usage and examples

### Yleaf: FASTQ (raw reads)
    
    python Yleaf.py -fastq raw_reads.fastq -out out -r 1 -q 20 -b 90 -t 16 -ref [hg19-hg38]
        
### Yleaf: BAM format
    
    python Yleaf.py -bam file.bam -out out -r 1 -q 20 -b 90 -t 1 -ref [hg19-hg38]

### Haplogroup prediction (w/output generated from Yleaf)

    python predict_haplogroup.py -input Output_files/ -out output.hg

3. See complete manual at the website:
    https://www.erasmusmc.nl/genetic_identification/resources/

4. Bug report

Please send an email at d.montielgonzalez@erasmusmc.nl if there is a problem getting the software up and running.

### References and Supporting Information
A. Ralf, et al., Yleaf: software for human Y-chromosomal haplogroup inference from next generation sequencing data (2018).

https://academic.oup.com/mbe/article/35/7/1820/4993044

