<img src="yleaf_logo.png" width="500" alt="Yleaf logo">

# Yleaf: software for human Y-chromosomal haplogroup inference from next generation sequencing data

#### Arwin Ralf, Diego Montiel Gonzalez, Kaiyin Zhong and Manfred Kayser

### Department of Genetic Identification 
#### Erasmus MC University Medical Centre Rotterdam, The Netherlands

## Requirements

    Operating system: Linux, macOS, or Windows (via standalone executable).
    Internet connection: when running for the first time for downloading the reference genome. Alternatively you 
                         can configure your own references.
    Data storage: For installation we recommend a storage capacity of > 8 GB. 

## Installation

### Option 1: Standalone executable (no Python required)

Download the pre-built binary for your platform from the [releases page](https://github.com/genid/Yleaf/releases):

- `yleaf-linux.tar.gz` — Linux x86-64
- `yleaf-macos.tar.gz` — macOS (Intel and Apple Silicon)
- `yleaf-windows.zip` — Windows x86-64

Extract and run the `yleaf` executable directly. samtools and bcftools are bundled — no external tools needed.

```bash
# Linux / macOS
tar xzf yleaf-linux.tar.gz
./yleaf/yleaf -h

# Windows (PowerShell)
Expand-Archive yleaf-windows.zip
.\yleaf\yleaf.exe -h
```

### Option 2: conda environment (recommended for source installs)

```bash
# first clone this repository to get the environment_yleaf.yaml
git clone https://github.com/genid/Yleaf.git
cd Yleaf
# create the conda environment — it will be called yleaf
conda env create --file environment_yleaf.yaml
# activate the environment
conda activate yleaf
# pip install the cloned yleaf into your environment. Using the -e flag allows you to modify the config file in your cloned folder
pip install -e .

# verify that Yleaf is installed correctly. You can call this command from any directory on your system
Yleaf -h 
```

### Option 3: manual install

```bash
# install python 3.7+ and libraries
apt-get install python3
pip3 install pandas numpy
# install external tools
sudo apt-get install minimap2 samtools bcftools
# clone the yleaf repository
git clone https://github.com/genid/Yleaf.git
cd Yleaf
pip install -e .

# verify that Yleaf is installed correctly
Yleaf -h 
```

After installation you can navigate to `yleaf/config.txt` and add custom paths for the reference genomes listed there. This prevents Yleaf from downloading them on first run and allows you to point to an existing reference. Positions are based on hg38, hg19, or T2T (hs1).

## Usage and examples

Here follow some minimal working examples of how to use Yleaf with different input files. There are additional options that can be used to tune how strict Yleaf is as well as options to get private mutations and a graph showing the positioning of predicted haplogroups in the haplogroup tree.

### BAM or CRAM format

    Yleaf -bam file.bam -o bam_output --reference_genome hg38
    Yleaf -cram file.cram -o cram_output --reference_genome hg38

For CRAM files the reference genome is usually auto-detected from the file header; `-rg` is only needed if auto-detection fails.

### FASTQ (raw reads)

    Yleaf -fastq raw_reads.fastq -o fastq_output --reference_genome hg38

### VCF input

    Yleaf -vcf variants.vcf.gz -o vcf_output --reference_genome hg38

### PLINK / SNP-array

    Yleaf -plink dataset.bed -o plink_output --reference_genome hg38

### With haplogroup tree visualisation and private mutations

    Yleaf -bam file.bam -o bam_output --reference_genome hg38 -dh -p

`-dh` generates a self-contained interactive HTML file with zoomable tree, per-haplogroup tabs, and PDF export.

### Ancient DNA samples

Use the `-aDNA` / `--ancient_DNA` flag when working with ancient DNA. This ignores G>A and C>T mutations, which are common post-mortem deamination artefacts and would otherwise be misinterpreted as derived alleles.

    Yleaf -bam ancient_sample.bam -o output --reference_genome hg38 --ancient_DNA

### Selecting a haplogroup tree

Yleaf supports multiple reference trees. Use the `--tree` flag to select one or more:

| Tree name   | Description                          |
|-------------|--------------------------------------|
| `yfull`     | YFull v14 (default)                  |
| `yfull_v10` | YFull v10.01 (legacy)                |
| `ftdna`     | FTDNA Y-haplotree                    |
| `isogg`     | ISOGG tree                           |

    # Single tree (default is yfull)
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
