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
# install python 3.7+ and libraries
apt-get install python3
pip3 install pandas
pip3 install numpy
# install minimap2 for FASTQ alignment
sudo apt-get install minimap2 
# install SAMtools
sudo apt-get install samtools
# clone the yleaf repository
git clone https://github.com/genid/Yleaf.git
# pip install the yleaf repository
cd Yleaf
pip install -e .

# verify that Yleaf is installed correctly. You can call this command from any directory on your system
Yleaf -h 
```
After installation you can navigate to the yleaf/config.txt folder and add custom paths for the files listed there. This will make sure that Yleaf does not download the files on the first go or downloads the files in the provided location. This allows you to use a custom reference if you want. Please keep in mind that custom reference files might cause other issues or give problems in combination with already existing data files. Positions are based on hg38, hg19 or T2T (hs1).

## Usage and examples
Here follow some minimal working examples of how to use Yleaf with different input files. There are additional options
that can be used to tune how strict Yleaf is as well as options to get private mutations as well as a graph showing 
the positioning of predicted haplogroups of all your samples in the haplogroup tree.

### Yleaf: FASTQ (raw reads)

    Yleaf -fastq raw_reads.fastq -o fastq_output --reference_genome hg38
        
### Yleaf: BAM or CRAM format

    Yleaf -bam file.bam -o bam_output --reference_genome hg19 
    Yleaf -cram file.cram -o cram_output --reference_genome hg38 

### With drawing predicted haplogroups in a tree and showing all private mutations

    Yleaf -bam file.bam -o bam_output --reference_genome hg19 -dh -p

### Ancient DNA samples

Use the `-aDNA` / `--ancient_DNA` flag when working with ancient DNA. This ignores G>A and C>T mutations, which are common post-mortem deamination artefacts and would otherwise be misinterpreted as derived alleles.

    Yleaf -bam ancient_sample.bam -o output --reference_genome hg38 --ancient_DNA

### Selecting a haplogroup tree

Yleaf supports multiple reference trees. Use the `--tree` flag to select one or more trees:

| Tree name   | Description                          |
|-------------|--------------------------------------|
| `yfull`     | YFull v14 (default)                  |
| `yfull_v10` | YFull v10.01 (legacy)                |
| `ftdna`     | FTDNA Y-haplotree                    |
| `isogg`     | ISOGG tree                           |

    # Single tree
    Yleaf -bam file.bam -o output --reference_genome hg38 --tree yfull

    # Multiple trees in one run (single pileup, per-tree prediction)
    Yleaf -bam file.bam -o output --reference_genome hg38 --tree yfull ftdna isogg

### Mixture analysis (forensic)

The `-mix` flag enables forensic mixture deconvolution: Yleaf identifies the contributing Y-haplogroups in a DNA mixture from multiple male donors.

    Yleaf -bam mixture.bam -o output --reference_genome hg38 --tree yfull -mix

Results are written to a `.mix` file per sample. Mixture analysis is tree-aware and supports all reference trees.

## Additional information

For a more comprehensive manual please have a look at the [yleaf_manual](yleaf_manual.pdf).

If you have a bug to report or a question about installation consider sending an email to 
 a.ralf at erasmusmc.nl or create an issue on GitHub.


### References and Supporting Information
A. Ralf, et al., Yleaf: software for human Y-chromosomal haplogroup inference from next generation sequencing data (2018).

https://academic.oup.com/mbe/article/35/5/1291/4922696
